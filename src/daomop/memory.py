import tracemalloc
import linecache

tracemalloc.start()

previous_snapshot = None


def display_top(key_type='lineno', limit=3, detail=False):
    snapshot = tracemalloc.take_snapshot()
    global previous_snapshot
    if previous_snapshot is not None:
        top_stats = snapshot.compare_to(previous_snapshot, key_type)
    else:
        top_stats = snapshot.statistics(key_type)

    if detail:
        print("Top %s lines" % limit)
        for index, stat in enumerate(top_stats[:limit], 1):
            frame = stat.traceback[0]
            # replace "/path/to/module/file.py" with "module/file.py"
            filename = os.sep.join(frame.filename.split(os.sep)[-2:])
            print("#%s: %s:%s: %.1f KiB"
                  % (index, filename, frame.lineno, stat.size / 1024))
            line = linecache.getline(frame.filename, frame.lineno).strip()
            if line:
                print('    %s' % line)

        other = top_stats[limit:]
        if other:
            size = sum(stat.size for stat in other)
            print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    msg = "TOTAL "
    print(f"{msg} allocated size: {total/1024**2:3.2f} MiB")
    previous_snapshot = snapshot
