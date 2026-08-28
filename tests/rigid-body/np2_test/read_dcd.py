import struct
import sys

def read_record(f):
    n_bytes = f.read(4)
    if len(n_bytes) < 4:
        return None
    n = struct.unpack('<i', n_bytes)[0]
    data = f.read(n)
    n2 = struct.unpack('<i', f.read(4))[0]
    assert n == n2, (n, n2)
    return data

path = sys.argv[1] if len(sys.argv) > 1 else "traj.dcd"

with open(path, 'rb') as f:
    header = read_record(f)
    magic = header[0:4]
    icntrl = struct.unpack('<20i', header[4:4+80])
    nframes = icntrl[0]
    has_cell = icntrl[19] if len(icntrl) > 19 else 0
    # some GENESIS dcd variants store cell flag differently; just try both
    title = read_record(f)
    natom_rec = read_record(f)
    natom = struct.unpack('<i', natom_rec[0:4])[0]
    print(f"magic={magic} nframes(icntrl0)={nframes} natom={natom}")

    frame = 0
    while True:
        pos_before = f.tell()
        rec = read_record(f)
        if rec is None:
            break
        # Detect whether this record is a unit-cell record (48 bytes = 6 doubles)
        if len(rec) == 48:
            cell = struct.unpack('<6d', rec)
            xrec = read_record(f)
            yrec = read_record(f)
            zrec = read_record(f)
        else:
            # no cell record; this record IS the X coords
            xrec = rec
            yrec = read_record(f)
            zrec = read_record(f)
            cell = None

        x = struct.unpack(f'<{natom}f', xrec)
        y = struct.unpack(f'<{natom}f', yrec)
        z = struct.unpack(f'<{natom}f', zrec)

        frame += 1
        if frame <= 3 or frame % 20 == 0 or True:
            print(f"--- frame {frame} ---")
            for i in range(natom):
                print(f"  atom {i+1}: {x[i]:9.3f} {y[i]:9.3f} {z[i]:9.3f}")
