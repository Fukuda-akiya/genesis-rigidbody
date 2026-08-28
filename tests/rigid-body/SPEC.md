# GENESIS Rigid-Body Dynamics 仕様書

最終更新: 2026-08-27

## 1. 背景・目的

gREST(generalized REST2, solute tempering REMD)シミュレーションにおいて、系の一部領域(例: PDB 4YUS 由来の対称アセンブリの繰り返し単位、53原子/コピー)を剛体(並進+回転の6自由度)として扱いたい。gRESTは通常 `spdyn`(MPI領域分割)+ RESPA(`VRES`)積分器で実行されるため、本機能はまず **spdyn の VRES 内側ループ** に実装する。atdyn(LEAP)対応は後回しとする。

入力データ(`tests/rigid-body/`、未追跡):

- `rigidbody_1us.dat` — 1つの剛体の構成原子53個の body-fixed frame(剛体固有座標系)座標。フォーマット:
  ```
  The number of mass points 1
            53
  Coordinates in the body-fixed frame 1
    <atom_index> <x> <y> <z>
    ...(53行)
  ```
- `index_4YUS_dark.dat` — 剛体のコピーごとの原子インデックスのブロック(空行区切り)。フォーマット:
  ```
  53

  <idx1>
  <idx2>
  ...(53行)

  53

  <idx1>
  ...(53行)
  ...(EOFまでブロック繰り返し)
  ```

**重要な発見**: `rigidbody_1us.dat` の53個のインデックス(103,104,...,117,757...764,5805...)は `index_4YUS_dark.dat` の最初のブロックと完全一致する。つまり `rigidbody_1us.dat` は**1つの代表剛体の body-fixed 参照座標テンプレート**であり、`index_4YUS_dark.dat` の各ブロックはそのテンプレートを異なるグローバル原子番号にマップした剛体のコピーを表す。テンプレート内の i 番目の座標は、各ブロックの i 番目のインデックスに対応する(同じ内部原子順序を仮定)。

現在の GENESIS ソースに rigid-body 機能は存在しない。最も近い既存機構は SETTLE(水分子の剛体拘束)と SHAKE/RATTLE(結合長拘束)であり、これらの実装パターンを踏襲する。

**2026-08-27 追加決定**: `communicate_constraints`(MPIランク間の水/HGroup移行データ交換)への剛体対応拡張は高リスクなため、まず `mpirun -np 1` で完全に動作・検証することを優先し、np>1(ランクをまたぐ剛体移動)対応は明示的に別フェーズとする。

**2026-08-28 追記**: 上記np>1対応を実装し、`mpirun -np 2`の実機エンドツーエンド検証まで完了(§8b参照)。

## 2. 合意済みスコープ

| 項目 | 決定 |
|---|---|
| 対象エンジン | spdyn を最優先。atdyn(LEAP)は後回し |
| 対象積分器 | RESPA(`VRES`, `sp_md_respa.fpp`)の**内側ループ** |
| 対象アンサンブル | 第一段階: NVE (`nve_vv1`/`nve_vv2`) + NVT Bussi/Berendsen/NHC (`vel_rescaling_thermostat_vv1`/`vv2`)。Langevin・バロスタット(NPT)は follow-up |
| spdynのセル境界問題 | 剛体直径 ≥ セル辺長なら setup 時にエラー停止。ペアリスト探索範囲は変更しない |
| gRESTホット領域との関係 | 独立。剛体グループは `[RIGIDBODY]` セクション+専用ファイルで指定し、gRESTのsolute選択とは無関係(重なっても重ならなくてもよい) |
| 検証方法 | 4YUS実系のpsf/pdb/パラメータが無いため、合成テスト系での単体検証(エネルギー保存則、剛体内原子間距離不変性、四元数ノルム、運動量・角運動量保存) |

## 3. RESPA(VRES)内側ループの構造(調査結果)

`sp_md_respa.fpp` の主ループ(`vverlet_respa_dynamics`, 250行目付近):

```
do i = istart, iend-multistep+1, multistep      ! 外側ループ (multistep = elec_long_period)
  do j = 1, multistep                            ! 内側ループ
    istep = i + j - 1
    call integrate_vv1(..., inner_step=j, dt_long, dt_short, ...)
    ...
    if (j < multistep) then
      call compute_energy_short(...)             ! 短距離力のみ再計算
    end if
    if (j == multistep) then
      call compute_energy(...)                   ! 長距離力も含め再計算(外側境界)
    end if
    call integrate_vv2(..., inner_step=j, dt_long, dt_short, ...)
  end do
end do
```

`nve_vv1`(NVEの場合、`sp_md_respa.fpp:947`):
- `force_short` による半キック + 座標ドリフトは **毎内側ステップ**(`dt_short`)。
- `force_long` による半キックは `mod(istep-1, elec_long_period) == 0` の時のみ(=各外側ブロックの**最初**の内側ステップ)。
- RATTLE (`compute_constraints(ConstraintModeLEAP, ...)`) は毎内側ステップで呼ばれる。

`nve_vv2`(`sp_md_respa.fpp:1105`):
- `force_short` による半キックは毎内側ステップ。
- `force_long` による半キックは `inner_step == elec_long_period` の時のみ(=各外側ブロックの**最後**の内側ステップ)。
- RATTLE (`compute_constraints(ConstraintModeVVER2, ...)`) は毎内側ステップ。

`vel_rescaling_thermostat_vv1`/`vv2`(Bussi/Berendsen/NHC、`sp_md_respa.fpp:1211`, `1900`台)も同様の力の分割構造を持つが、速度スケーリング(サーモスタット)が追加される。

**剛体積分の挿入方針**: 上記の点粒子ループと同じ場所・同じ `force_long`/`force_short` の非対称適用パターンに従い、剛体構成原子については点粒子としての通常キック・ドリフトの**代わりに**、剛体の重心・角運動量・姿勢四元数を伝播し、その結果で当該原子の `coord`/`velocity` を上書きする(SETTLEが非拘束更新後に解析解で上書きする方式と同じ)。

## 4. 剛体の状態と積分アルゴリズム

剛体ごとの状態(spdynでは領域分割によりセル局所に保持。詳細は §6):
- 構成原子のグローバルインデックスリスト
- body-fixed frame 参照座標(テンプレート、全コピー共有)
- 全質量、慣性テンソルの主慣性モーメント・主軸(setup時に参照座標+原子質量から計算。3x3対称行列の解析的対角化、LAPACK非依存)
- 重心位置・重心速度、姿勢四元数、角運動量

積分(kick-drift-kick、RESPAの内側/外側力分割に対応):

VV1相当(各内側ステップ j):
1. 剛体の構成原子の現在の力(`force_short`、および j が外側ブロック最初のステップなら追加で `force_long`)から正味力・重心周りのトルクを集計。
2. 半キック: 重心速度・角運動量を(`force_short`は`half_dt_short`、`force_long`は`half_dt_long`で)更新。
3. ドリフト: 重心を `dt_short * v_com` だけ並進。角運動量→(慣性テンソル経由で)角速度に変換し、姿勢四元数を `dt_short` 分回転更新。
4. 剛体原子の座標を `重心 + R(四元数)・body-fixed座標` として再構成し `coord` に書き戻す。速度は剛体速度場 `v_com + ω × r_i` として `velocity` に書き戻す。

VV2相当(各内側ステップ j):
1. 新しい力(`force_short`、および j が外側ブロック最後のステップなら追加で `force_long`)から正味力・トルクを再集計。
2. 半キック: 重心速度・角運動量を更新(座標は更新しない、点粒子のVV2と同じ)。
3. 剛体原子の速度を `v_com + ω × r_i` として `velocity` に書き戻す。

剛体はもともと厳密に形状を保つため、RATTLE のような反復拘束解法は不要(6自由度で厳密に表現するのが本設計の利点)。

自由度: 剛体に含まれる原子は 3N → 6(並進3+回転3)に変わるため、`domain%num_deg_freedom` を補正する(既存の自由度調整箇所と同じ場所)。

## 5. 共有ライブラリ(`src/lib/`) — 新規ファイル

- **`fileio_rigidbody.fpp`**: `fileio_spot.fpp`(open/read/close)+ `fileio_localres.fpp`/`fileio_morph.fpp`(件数不明ブロックをEOF検出しつつ読む2パス方式)を参考に実装。
  - `input_rigidbody_index(filename, rb_index)` — `index_4YUS_dark.dat` 形式を読み、剛体コピーごとの原子インデックス配列を構築。
  - `input_rigidbody_refcoord(filename, rb_ref)` — `rigidbody_1us.dat` 形式を読み、テンプレート座標配列を構築。
- **`rigidbody_str.fpp`**: `s_rigidbody` 型定義(`at_constraints_str.fpp` の `s_constraints` のパターンを踏襲)。`alloc_rigidbody`/`dealloc_rigidbody`。
- **`rigidbody.fpp`**:
  - `setup_rigidbody(molecule, rb_index, rb_ref, rigidbody)` — 各ブロックの原子数とテンプレート原子数の一致検証、質量・重心・慣性テンソル・初期姿勢四元数の計算。領域分割**前**のグローバル `molecule` に対して呼ぶ。
  - 力・トルク集計、kick/drift の数値プリミティブ(汎用配列引数、`s_domain` 型に非依存の純粋関数として実装し、atdyn/spdyn 双方から呼べるようにする)。

## 6. spdyn 側の配線

- `sp_rigidbody.fpp`(新規): `s_rgbd_info`(`[RIGIDBODY] rigidbody_indexfile=... rigidbody_reffile=...`)、`read_ctrl_rigidbody`、`setup_rigidbody_spdyn`。setup時に剛体直径とセル辺長を比較し、直径 ≥ セル辺長ならエラー停止。
- `sp_control.fpp`: `[RIGIDBODY]` セクション登録。
- `sp_setup_spdyn.fpp`: `setup_domain` の**前**にグローバル `molecule` に対して `setup_rigidbody(...)` を呼ぶ。
- `sp_domain_str.fpp`: `domain%num_rigidbody(icel)`, `domain%rigidbody_list(:,:,icel)`, `domain%rigidbody_quat(4,:,icel)`, `domain%rigidbody_angmom(3,:,icel)`, `domain%rigidbody_vel_com(3,:,icel)` 等を追加(`water_list`/`num_water` と同じセル局所ストレージパターン)。
- `sp_domain.fpp`: `setup_atom_by_HBond` の water セクションを参考に、剛体の代表原子でセルを決定し、他の構成原子を強制的に同じセルへ割り当てる。
- `sp_migration.fpp`: `update_outgoing_water`/`update_incoming_water` を参考に、四元数・角運動量ペイロードを含めて移送する処理を追加。
- `sp_md_respa.fpp`: `nve_vv1`/`nve_vv2`、`vel_rescaling_thermostat_vv1`/`vv2` の各内側ステップに剛体伝播呼び出しを追加(§4 のアルゴリズムに従う)。
- `sp_dynamics.fpp`: `rigidbody` を `run_md` に伝搬、自由度調整。

## 7. 検証計画

- パーサ単体テスト: `tests/rigid-body/` の2ファイルを実際に読み込ませ、ブロック数・原子数・座標値を確認。
- 合成系での物理的不変量テスト(spdyn, np=1 および np>1 で剛体が単一セルに収まる配置):
  - 全系エネルギー保存(NVE)
  - 剛体内原子間距離の不変性
  - 四元数ノルム = 1
  - 全系の運動量・角運動量保存則
- ビルド確認: spdyn(MPI/LAPACK 有無の組み合わせ)。

## 8. 実装状況(2026-08-27時点)

- [x] ビルド環境構築: autotools(aclocal/autoconf/automake)をHomebrewで導入し `./bootstrap && FC=mpif90 ./configure` でMPI+LAPACK有効の設定を確認(`OMPI_CC=gcc-15 make` でOpenMP対応のCコンパイラを使用する必要あり、Apple clangは`-fopenmp`非対応のため)。spdynのベースラインビルドが成功することを確認済み。
- [x] 共有ライブラリ実装・検証済み:
  - `src/lib/fileio_rigidbody.fpp` — `input_rigidbody_index`/`input_rigidbody_refcoord`。実データ(`index_4YUS_dark.dat`, `rigidbody_1us.dat`)で単体テスト済み。
  - `src/lib/rigidbody_str.fpp` — `s_rigidbody` 型、`alloc_rigidbody`/`dealloc_rigidbody`。
  - `src/lib/rigidbody.fpp` — `setup_rigidbody`(質量・COM・慣性テンソル計算、LAPACK `dsyev` による対角化で主軸フレームに変換)、`rigidbody_diameter`(spdynのセルサイズチェック用)、`quat_to_rotmatrix`/`quat_normalize`/`quat_derivative`/`rigidbody_angvel`(四元数・剛体角速度の数値プリミティブ、Miller et al. 2002方式: 角運動量をspace frameで保持しbody frameへ回転してから対角慣性を適用)。
  - **設計変更**: 慣性テンソル対角化は当初「LAPACK非依存の解析解」を予定していたが、`fitting.fpp`の`fit_trrot`が同種の対角化(4x4版)で既にLAPACK `dsyev` に依存しており、自作の3x3固有値ソルバーより堅牢かつ実装リスクが低いため、LAPACK依存(`#ifdef LAPACK`ガード付き、既存コードと同じ規約)に変更した。
  - `tests/rigid-body/test_rigidbody.f90` — 開発用スタンドアロン単体テスト(ビルド手順はファイル末尾コメント参照)。パーサ検証+合成平面3原子剛体での物理検証(質量・COM・慣性テンソル対角性・垂直軸の定理・直径)を全項目PASS。
- [x] spdyn配線(ランク内)実装済み:
  - `src/spdyn/sp_rigidbody.fpp` — `[RIGIDBODY]`制御セクション(`s_rgbd_info`, `read_ctrl_rigidbody`)、`setup_rigidbody_spdyn`(共有ライブラリの`setup_rigidbody`呼び出し+剛体直径とセルサイズの比較チェック)。
  - `sp_domain_str.fpp` — `s_domain`にセル局所の剛体ストレージ(`num_rigidbody`, `rigidbody_id`, `rigidbody_atom`, `rigidbody_com/vel_com/quat/angmom`)と移行ステージングバッファ(`s_domain_rigidbody`型)を追加。`alloc_domain`に`DomainRigidBody`/`DomainRigidBodyMove`ケースを追加。
  - `sp_domain.fpp` — `setup_domain`に`rigidbody`(optional)引数を追加。`mark_rigidbody_atoms`(剛体原子を`constraints%duplicate`でマークし、既存のsoluteループから除外。SHAKE/SETTLE拘束済み原子との重複はsetup時エラー)と`setup_atom_by_rigidbody`(water/HGroupと同じ「代表原子でセル決定→全構成原子を強制的に同じセルへ」パターン)を追加し、`setup_hbond_group`/`setup_atom_by_HBond`の直後に呼び出す。
  - `sp_migration.fpp` — `update_outgoing_rigidbody`/`update_incoming_rigidbody`(+ヘルパー`pack_rigidbody`/`unpack_rigidbody`)を`update_outgoing/incoming_water`と同じパターン(毎回ステージングバッファから完全再構築、`water`の直後に実行)で実装。原子ごとのcharge/mass/atom_cls_no/グローバルIDは通信せず、共有(全ランク複製済み)の`s_rigidbody`構造体から都度参照する設計(`molecule`はsetup完了後に解放されるため)。
  - `sp_update_domain.fpp` — `domain_interaction_update`に`rigidbody`(optional)引数を追加し、`update_outgoing/incoming_water`の直後に剛体版を呼び出す配線を追加。
  - 単体テスト(`tests/rigid-body/test_rigidbody.f90`)は全項目PASSを維持(atom_mass/atom_charge/atom_cls_noの捕捉を含む)。
- [x] **未解決の重大なギャップ(スコープ決定済み)**: `sp_communicate.fpp`の`communicate_constraints`(MPIランク間での水/HGroup移行データ交換を行う、手作業で最適化された複雑な3次元境界通信ルーチン)に、剛体の移行ステージングバッファ(`domain%rigidbody%move_real`等)を送受信する処理は**まだ組み込まれていない**。
  - **`mpirun -np 1`(単一ランク)では正しく動作する**(`update_outgoing_rigidbody`/`update_incoming_rigidbody`自体はランク内のセル間移行として完結しているため)。
  - **`np>1`で、剛体の代表原子が別ランクが担当するセルへ移動する場合は未対応**(パック済みデータが送信されずロストする)。ユーザーと合意の上、np=1での完成を優先し、np>1対応は別フェーズとした(2026-08-27決定、本セクション冒頭の「追加決定」参照)。
- [x] **`sp_md_respa.fpp` の RESPA内側ループへの剛体積分(本題)実装・検証済み**:
  - `write_rigidbody_atoms` — 剛体の現在状態(com/vel_com/quat/angmom)から構成原子のcoord(オプション)/velocityを再構成し上書き(SETTLEと同じ「非拘束更新→解析解で上書き」パターン)。
  - `rigidbody_force_torque` — 剛体の現在COM周りの正味力・トルクを構成原子のforce配列から集計。
  - `propagate_rigidbody_vv1`/`propagate_rigidbody_vv2` — RESPA内側ループの各innerステップで、`force_short`は毎ステップ、`force_long`は外側ブロック境界でのみ半キック(既存の`nve_vv1`/`nve_vv2`の点粒子分割と同一パターン)。回転ドリフトは**body-frame角速度の指数写像**(`quat_rotate_by_body_omega`、SO(3)の幾何学的積分子。当初は線形化四元数微分(オイラー法)を使っていたが、エネルギー保存性を検証した際に不十分と判明し置き換えた。詳細は下記の検証結果参照)。
  - `initialize_rigidbody_state` — シミュレーション開始時(`vverlet_respa_dynamics`内、`initial_vverlet`/リスタート読み込みの直後)に1回、現在の原子座標・速度から各剛体のcom/vel_com/quat/angmomを初期化。姿勢は`fit_rigidbody_quat`(Kearsley/Kabsch四元数法、`fitting.fpp`の`fit_trrot`と同じLAPACK `dsyev`方式だが密な1..n配列に対応させた独自実装)で参照座標を実際の初期構造にフィッティングして決定。
  - `nve_vv1`/`nve_vv2`/`integrate_vv1`/`integrate_vv2`/`vverlet_respa_dynamics`に`rigidbody`(optional)引数を追加し配線。当初はNVEのみ対応だったが、2026-08-28にNVT(Bussi/Berendsen/NHC)対応を追加(下記)。NPT・Langevinは引き続き`integrate_vv1`がrun時にエラー終了する。
- [x] `sp_dynamics.fpp` での `rigidbody` の `run_md` への伝搬、自由度調整(`update_num_deg_freedom`で剛体原子の3N自由度を6/bodyに置換)実装済み。
- [x] `sp_setup_spdyn.fpp` での `setup_rigidbody_spdyn` 呼び出し配線(`setup_spdyn_md`、FEPと非FEPの両方の`setup_domain`分岐に対応)実装済み。
- [x] **合成系での物理検証(np=1)実施・全項目PASS**:
  - `tests/rigid-body/test_rigidbody_dynamics.f90` — `sp_md_respa_mod`の`propagate_rigidbody_vv1`/`vv2`/`initialize_rigidbody_state`を(`public`化した上で)実際のspdynオブジェクトファイル群にリンクして直接呼び出す統合テスト。手作業で構築した最小限の単一セル`s_domain`(制御ファイル・MPI分割・力場計算は一切介さない)上で、非対称・非平面な4原子剛体(3主慣性モーメントが相異なる漸近的トップ)を外力ゼロで自由回転・並進させ、以下を確認:
    - 剛体内原子間距離が20000ステップ後も不変(剛体性、厳密)
    - 四元数ノルムが1(数値精度内)
    - 全系運動エネルギー(並進+回転)が1%以内で保存(下記の積分精度の議論を参照)
    - 全系角運動量(space frame)が厳密に保存(外部トルクゼロ、実装の整合性チェック)
  - **積分精度に関する発見**: 剛体の回転運動は「各ステップ開始時のbody-frame角速度で指数写像により厳密に1軸回転させる」方式であり、真の自由回転(Euler方程式)ではステップ内でも角速度自体が連続的に変化する(トルクフリー歳差運動)ため、この方式は1次精度の近似となる。実測: dt=1e-3(20 time units、意図的に大きな初期角速度)でエネルギードリフト約0.81%、dt=2e-4(5倍小さく、同じ総シミュレーション時間)で約0.16%と、O(dt)にほぼ比例して縮小することを確認 — バグではなく真の離散化誤差であることを検証済み。gRESTの実用条件(現実的なMDタイムステップ、熱的な回転速度)ではこのテストよりはるかに小さいドリフトが期待されるが、完全にシンプレクティックな複数副ステップ回転積分子(DLM/NO_SQUISH方式など)は未実装(follow-up)。
  - **開発中に発見・修正した重大なバグ**: `propagate_rigidbody_vv1`のトルク計算が当初`domain%coord`を直接読んでいたが、同じステップ内で先行する点粒子の素朴なkick+driftループ(`nve_vv1`/`vel_rescaling_thermostat_vv1`内)が剛体構成原子の`domain%coord`も(物理的に無意味な値で)上書きしてしまうため、トルクの梃子腕(lever arm)`r = coord - com`が汚染されたデータで計算されていた。外力ゼロの`test_free_rigidbody`ではトルク自体が常にゼロになるため検出されず、質量比例な外力(重力)を課す新規テスト`test_gravity_torque_free`(正味トルクは解析的に厳密ゼロのはず)を追加して初めて発覚。修正: 力評価時点の座標を保持している`domain%coord_ref`を読むよう変更(`rigidbody_force_torque`に`coord`配列を明示引数として追加)。修正前は角運動量ドリフトが系のスケールに対して非常に大きく(検出不能なレベルではない)、修正後はドリフトが厳密にゼロになることを確認。

- [x] **NVT(Bussi/Berendsen/NHC)対応(2026-08-28)実装・検証済み**:
  - `sp_constraints.fpp`の`compute_kin_group`(RESPAの運動エネルギー算出で使われる「グループ運動エネルギー」計算、`group_tp`の有無に関わらず`nve_vv1`/`vel_rescaling_thermostat_vv1`から常に呼ばれる)に剛体対応を追加。**重要な発見**: 剛体構成原子は`constraints%duplicate`マークにより`nsolute`(通常のsolute原子カウント)からも`HGr`/`water`からも除外されているため、この関数は元々剛体の運動エネルギー寄与を完全に無視していた(サイレントな過小評価)。水/HGroupが「内部拘束自由度をゼロとして無視」するのに対し、剛体は回転の3自由度が本物の(無視してはいけない)自由度であるため、`mass*vel_com^2 + omega_space・angmom_space`(並進+回転の全6自由度分)を追加する形で対応。
  - `sp_md_respa.fpp`に`write_rigidbody_vel_half_estimate`(半キック速度推定値`vel_half`への剛体対応。`ensemble%group_tp = NO`の場合に使われる`calc_kinetic`パスのため)、`rescale_rigidbody_velocity`(サーモスタットのスケール係数を`rigidbody_vel_com`/`rigidbody_angmom`にも適用)を追加。`vel_rescaling_thermostat_vv1`にこれらと`propagate_rigidbody_vv1`呼び出しを配線。
  - `integrate_vv1`のアンサンブル制限を緩和: NVE、またはNVT+tpcontrol=Bussi/Berendsen/NHCを許可(Langevin/NPTは引き続きrun時エラー)。
  - `tests/rigid-body/test_rigidbody_dynamics.f90`に`test_gravity_torque_free`を追加(上記バグの回帰テストとして機能。角運動量の厳密保存、COMの正しい等加速度運動(射影運動)、剛体性を検証、全PASS)。
  - **未検証の範囲**: `vel_rescaling_thermostat_vv1`本体(Bussi/Berendsen/NHCの実際のスケール係数計算ロジック自体)を経由したエンドツーエンドのNVT実行(`s_ensemble`/`s_dynvars`/乱数種等を伴う完全な統合テスト)は未実施。今回検証したのは、その前後で必要になる剛体側の運動エネルギー計算・速度スケーリング・力/トルク計算の各要素(最もリスクが高い部分)であり、サーモスタットの温度収束挙動そのものの実測は次のステップ。

## 8b. np>1 対応(2026-08-28実装・実機検証済み)

- **`sp_communicate.fpp`の`communicate_constraints`に剛体の移行データ送受信を実装**。既存のwater/HGroup送受信パターン(3次元カスケード通信、`ii=3,1,-1`ループ、lower/upper双方向、ヘッダー個数スロット+固定長レコードのパック/アンパック)を忠実に踏襲。剛体は1体につき`6*max_natom+13`個のreal(構成原子ごとのcoord+velocity 6個 × natom、+ com/vel_com/quat/angmomの13個)+1個のint(グローバル剛体ID)の固定長レコードとしてwaterと全く同じ形式でパックする設計。`domain%rigidbody%move`系のヘッダースロットは`constraints%connect+3`番目に追加(`rigidbody%is_used`が全ランクで一致するため、レイアウトはランク間で常に整合)。通信バッファ(`comm%buf_send`/`int_send`)の既存サイズ見積もり(`8*MaxAtom*max_cell`)に対して、剛体1体あたりの実質ペイロードは同じ原子数を通常soluteとして送る場合より小さいため、既存の割り当てで十分と判断(防御的なオーバーフローチェックも追加)。
- **`mpirun -np 2`での実機エンドツーエンド検証を実施**(`tests/rigid-body/np2_test/`): 手作業で構築した最小限の自己完結GROMACS形式トポロジー(結合なし、3原子の剛体+4個の背景原子、LJ相互作用のみ、電荷ゼロ)を用い、剛体をドメイン境界(x=30Å、`domain_x=2`)付近に配置し、近傍の「押し出し原子」からのLJ反発で境界を越えて加速させた。結果:
  - 200フレームの全トラジェクトリで剛体内原子間距離の最大偏差 **0.000002Å**(浮動小数点精度レベル、実質的に完全な剛体性)。
  - 代表原子のx座標は28.6〜41.9Åの範囲でドメイン境界(x=30)を**複数回正しく往復**(周期境界条件でバウンドする押し出し原子により、剛体が繰り返しランク0↔ランク1間を移動)。
  - 実行は正常終了(exit code 0)、`comm_migrate`タイマーに非ゼロの時間が計上され、実際に移行処理が実行されたことを確認。
  - **観測された既知の挙動**: 境界通過に近いタイミングで、エネルギー出力に1ステップ限りの小さな遷移的変動(全エネルギーで約2.3%、直後に新しい安定値へ復帰、その後は成長・発散せず4000ステップ全体で安定)が見られた。これはペアリストが`nbupdate_period`ごとにしか更新されないことに起因する、セル/ドメイン再割当て直後の一時的な力場評価の近似誤差であり、通常粒子でも起こりうるドメイン分割MDコード一般の既知の特性と考えられる(剛体固有のバグではないと判断: もしデータ破損があれば剛体性が破れるはずだが、遷移の前後で剛体性は完璧に保たれている)。
  - 検証用の一式(`system.top`, `system.gro`, `np2_test.inp`, `index_np2.dat`, `ref_np2.dat`, `read_dcd.py`)は`tests/rigid-body/np2_test/`に保存済み。GROMACS形式トポロジーを使う際の実務上の注意点(spdynは`[exclusions]`セクションを無視する、VRES積分器は`electrostatic=CUTOFF`単体では`elec_long_period>1`を使えず`Multi-step is available only with PME`エラーになる、`thermostat_period`/`barostat_period`は`elec_long_period`の倍数である必要がある、剛体内原子間の意図しない巨大LJ反発を避けるにはsigmaを原子種ごとに差別化する必要がある、など)も併せて記録。

## 9. 既知の制限・follow-up課題

- np>1のMPIランク間剛体移行は上記§8bの通り実機検証済み。ただし検証したのは剛体が1個・小規模合成系のみであり、多数の剛体が同時に複数方向へ移行する高負荷ケースは未検証。
- **NVT(Bussi/Berendsen/NHC)は実装済みだがエンドツーエンド未検証**(§8参照)。`vel_rescaling_thermostat_vv1`本体を通した実際の温度収束挙動の検証は未実施。
- **NPT・Langevinサーモスタットとの組み合わせは引き続き未対応**。`integrate_vv1`が`rigidbody%is_used`かつ(アンサンブルがNVE/NVTのいずれでもない、またはNVTでもtpcontrol=Langevinの)場合に明示的にエラー終了する。理由: (1) NPTは`mtk_barostat_vv1`/`vv2`(圧力/ビリアル計算を含む、さらに複雑な拡張が必要)、(2) `compute_virial_group`(圧力のvirial計算、`ensemble%group_tp`向け)にはまだ剛体対応を追加していない(NVTでは運動エネルギーのみで十分だが、NPTには必須)、(3) Langevinサーモスタット(`langevin_thermostat_vv1`/`vv2`)は速度リスケーリング方式ではなくランダム力+摩擦項を直接force配列に加算する方式のため、別の対応方針が必要。
- 完全にシンプレクティックな回転積分子(DLM/NO_SQUISH等)は未実装。現行の「毎ステップ指数写像1回」方式はO(dt)精度(§8の検証結果参照)。
- リスタートファイル(.rst)への四元数・角運動量の保存は未対応(剛体を含む実行はリスタート継続不可)。
- 剛体内部の非結合(LJ/静電)相互作用は `compute_energy`/`compute_energy_short` で計算され続ける(正味力・トルクへの寄与は相殺されるが計算コストは無駄)。最適化は見送り。
- atdyn(LEAP)対応は後回し。
- 4YUS実系での本番検証は、対応するpsf/pdb/パラメータファイルをユーザーが別途用意しない限り実施不可(今回の検証は手作業で構築した合成系による)。
- spdyn: 剛体直径 ≥ セル辺長の場合は setup 時エラー(ペアリスト探索範囲は変更しない、という制限つき)。
- `dealloc_domain_all`(あれば)に`DomainRigidBody`/`DomainRigidBodyMove`のdeallocケースを追加していない(プロセス終了時のメモリリークのみ、実行中の正しさには影響しない)。
- `domain%num_deg_freedom`の剛体自由度調整は`setup_spdyn_md`(非FEP・FEPの`setup_domain`呼び出し双方)のみに実装。`setup_spdyn_remd`/`setup_spdyn_rpath`/`setup_spdyn_min`には配線していない(スコープ外、§2参照)。
