#!/usr/bin/env python3
"""
This script is the entrypoint to collecting and summarizing run data for
our pipeline. It is focused primarily on timing information and cluster
utilization.
"""

# import sys
import io
import os
import logging
import math
import glob
import csv
import re
from datetime import datetime, timedelta

import typer
from typing import Optional, Iterable, List
from dataclasses import dataclass

log = logging.getLogger(__name__)

app = typer.Typer()

# TODO: Make these better
seq_roots = ["/net/seq/data2/sequencers", "/net/seq/data/sequencers"]
flowcell_roots = ["/net/seq/data2/flowcells", "/net/seq/data/flowcells"]
agg_roots = ["/net/seq/data2/aggregations/", "/net/seq/data/aggregations/"]


@dataclass
class Event():
    """
    An Event represents a step or task.
    It must have a name and a start_time - duration and realtime may be 0 if
    there is no logical duration/end-time for the event.
    """
    source: Optional[str]
    name: str
    start_time: datetime
    duration: timedelta
    realtime: timedelta
    cpus: Optional[float]
    memory: Optional[int]
    bytes_read: Optional[int]
    bytes_written: Optional[int]

    @staticmethod
    def from_file(event_name: str, filename: str) -> "Event":
        """
        Create an Event based on the file's mtime
        Duration and realtime will be 0.
        """
        return Event(
            source=filename,
            name=event_name,
            start_time=datetime.fromtimestamp(os.path.getmtime(filename)),
            duration=timedelta(0),
            realtime=timedelta(0),
            cpus=None,
            memory=None,
            bytes_read=None,
            bytes_written=None,
        )

    @staticmethod
    def from_file_pair(event_name: str,
                       start_filename: str,
                       end_filename: str) -> "Event":
        """
        Create an Event based on two files' mtime
        """
        start_time = datetime.fromtimestamp(os.path.getmtime(start_filename))
        end_time = datetime.fromtimestamp(os.path.getmtime(end_filename))
        return Event(
            name=event_name,
            start_time=start_time,
            duration=end_time - start_time,
            realtime=end_time - start_time,
            cpus=None,
            memory=None,
            bytes_read=None,
            bytes_written=None,
            source=end_filename,
        )

    @property
    def end_time(self) -> datetime:
        """Calculates the time that this event ended"""
        return self.start_time + self.duration

    def __str__(self) -> str:
        return self.to_readable()

    def to_readable(self) -> str:
        return "\t".join([self.name, str(self.duration), str(self.realtime)])


@dataclass
class Events():
    """
    Basically just a list of Event, plus some utility methods.
    """
    events: List[Event]

    def start_time(self) -> datetime:
        """ Gets start of the first event """
        return min(e.start_time for e in self.events)

    def end_time(self) -> datetime:
        """ Get the end of the last event """
        return max(e.end_time for e in self.events)

    def elapsed_time(self) -> timedelta:
        """
        Computes the time between the first start and the last end
        """
        return self.end_time() - self.start_time()

    def merged(self, other: 'Events') -> 'Events':
        new_list = self.events + other.events
        return Events(events=new_list)

    def to_tsv(self, print_header=True) -> str:
        """ Render the event list as a TSV """
        def format_duration(dur):
            return round(dur / timedelta(hours=1), 2)

        def format_cpus(cpus):
            return round(cpus, 2) if cpus else None

        output = io.StringIO()
        writer = csv.DictWriter(
            output,
            fieldnames=[
                "source",
                "name",
                "start_time",
                "end_time",
                "duration_hours",
                "realtime_hours",
                "memory_bytes",
                "cpus",
                "bytes_read",
                "bytes_written",
            ],
            delimiter="\t",
        )
        if print_header:
            writer.writeheader()
        for event in self.events:
            writer.writerow({
                "source": event.source,
                "name": event.name,
                "duration_hours": format_duration(event.duration),
                "start_time": event.start_time,
                "end_time": event.end_time,
                "realtime_hours": format_duration(event.realtime),
                "memory_bytes": event.memory,
                "cpus": format_cpus(event.cpus),
                "bytes_read": event.bytes_read,
                "bytes_written": event.bytes_written,

            })
        return output.getvalue()


def parse_trace_datetime(time: str) -> datetime:
    """Parses the datetime used by trace files"""
    # ISO format datetime
    return datetime.strptime(
        time, r"%Y-%m-%d %H:%M:%S.%f")


# A time regex to parse times like '3h 2m' or '12m 18s'
time_re = re.compile(
    r'(?:(?P<hours>\d+)h)?'
    r'\s*'
    r'(?:(?P<minutes>\d+)m)?'
    r'\s* (?:(?P<seconds>\d+)s)?',
    flags=re.VERBOSE)


def parse_trace_duration(time: str) -> timedelta:
    """
    Parses a duration in trace.txt format into a timedelta
    Raises an exception if cannot parse
    """
    match = time_re.match(time)
    if not match:
        raise Exception(f"Could not parse timedelta '{time}'")
    groups = match.groupdict()
    hours = int(groups['hours'] or 0)
    minutes = int(groups['minutes'] or 0)
    seconds = int(groups['seconds'] or 0)
    return timedelta(hours=hours,
                     minutes=minutes,
                     seconds=seconds,
                     )


def try_parse_bytes(size: str) -> Optional[int]:
    try:
        return parse_bytes(size)
    except ValueError:
        return None


def parse_bytes(size: str) -> int:
    """
    Returns the number of bytes
    Raises if cannot parse
    """
    power_table = {
        "B": 0,
        "KB": 1,
        "MB": 2,
        "GB": 3,
        "TB": 4,
    }
    parts = size.strip().upper().split()
    if len(parts) < 1:
        raise Exception("Empty string is not a bytes spec")
    if len(parts) > 2:
        raise Exception(f"bytes value looks malformed: {size}")
    num = float(parts[0])
    if len(parts) == 1:
        return math.floor(num)

    base = parts[1]
    if base not in power_table:
        raise Exception(f"Unknown suffix: {base}")
    return math.floor(num * (1024 ** power_table[base]))


def try_parse_percent(val: str) -> Optional[int]:
    try:
        return parse_percent(val)
    except ValueError:
        return None

def parse_percent(val: str) -> float:
    return float(val.replace("%","").strip()) / 100


def parse_trace_file(filename: str, source=None) -> List[Event]:
    if source is None:
        source = filename
    with open(filename, 'r') as f:
        return parse_trace(f.read(), source=source)


def parse_trace(trace: str, source=None) -> List[Event]:
    """
    Parses the contents of a trace file into a list of Events
    """
    return [
        Event(source=source,
              start_time=parse_trace_datetime(row['submit']),
              name=row['name'],
              duration=parse_trace_duration(row['duration']),
              realtime=parse_trace_duration(row['realtime']),
              memory=try_parse_bytes(row.get('rss') or row['peak_rss']),
              cpus=try_parse_percent(row[r'%cpu']),
              bytes_read=try_parse_bytes(row['rchar']),
              bytes_written=try_parse_bytes(row['wchar']),
              )
        for row in csv.DictReader(trace.splitlines(), delimiter="\t")
    ]


def parse_run_file(run_script: str) -> List[str]:
    """
    Returns the directories from a run_*.bash file (the kind created by
    alignprocess.py)
    """
    out = []
    with open(run_script, 'r') as f:
        for line in f.readlines():
            words = line.split()
            if len(words) >= 2 and words[0] == "cd":
                out.append(words[1])
    return out

def get_elapsed_time(directory: str) -> timedelta:
    """
    Returns the difference between the newest and oldest items in a directory
    Does not recurse into sub-directories.
    """
    paths = os.listdir(path=directory)
    path_times = [
        timedelta(seconds=os.path.getmtime(os.path.join(directory, f)))
        for f in paths]
    return max(path_times) - min(path_times)


def get_trace_contents(directory: str) -> str:
    """ Gets the content of the most recent trace in the directory """
    trace_filename = os.path.join(directory, "trace.txt")
    with open(trace_filename, mode="r", encoding="utf8") as f:
        return f.read()


def get_seq_events(seq_dir: str) -> Events:
    events: List[Event] = []
    seq_events = [
        ("Sequencing", "Recipe", "RTAComplete.txt"),
        ("Setup", "processing.json", "run_bcl2fastq_2.sh"),
        ("BarcodeCount", "bc-*.o*", "bc-*.e*"),
        ("bcl2fastq", "u-*.o*", "u-*.e*"),
        ("Demux", "queuedemux-*.o*", "c-*.e*"),
        ("Copy", "c-*.o*", "c-*.o*"),
        ("Collate", "c-*.o*", "collate-*.o*"),
    ]
    for (name, start_file, end_file) in seq_events:
        events.append(Event.from_file_pair(
            name,
            glob_one(seq_dir, start_file),
            glob_one(seq_dir, end_file),
        ))

    return Events(events=events)


def get_flowcell_events(fc_dir: str) -> Events:
    events: List[Event] = []

    def first_mtime(*glb_parts: str) -> datetime:
        glb = os.path.join(fc_dir, *glb_parts)
        return datetime.fromtimestamp(os.path.getmtime(
            min(glob.glob(glb),
                key=os.path.getmtime)
        ))

    def last_mtime(*glb_parts: str) -> datetime:
        return datetime.fromtimestamp(os.path.getmtime(
            max(glob.glob(os.path.join(fc_dir, *glb_parts)),
                key=os.path.getmtime)
        ))
    # TODO: Where are the log files?
    align_start = first_mtime("Project_*", "Sample_*", "align*", "*run.bash")
    align_end = last_mtime("Project_*", "Sample_*", "align*", "trace.txt")
    align_time = align_end - align_start

    events.append(
        Event(name="Alignment", start_time=align_start,
              duration=align_time, realtime=align_time)
    )
    agg_dirs = parse_run_file(os.path.join(fc_dir, "run_aggregations.bash"))

    last_agg_trace = max(
        [os.path.join(d, "trace.txt")
         for d in agg_dirs
         if os.path.exists(os.path.join(d, "trace.txt"))],
        key=os.path.getmtime
    )
    events.append(
        Event.from_file_pair("Aggregation",
                             os.path.join(fc_dir, "run_aggregations.bash"),
                             last_agg_trace)
    )

    return Events(events)

# def get_agg_dirs(agg_ids: Iterable[int]) -> Iterable[str]:
#     # Doesn't work
#     # brace expansion {123,124,125} doesn't work with glob
#     needle = (
#         "aggregation-{"
#         + ",".join(str(i) for i in agg_ids)
#         + "}"
#     )
#     return glob.glob(os.path.join(agg_root, "LN*", needle))


def glob_one(*parts: str) -> Optional[str]:
    """
    A lot of times, we expect to find exactly one file matching a glob.
    This makes those scenarios faster.
    If two or more exist, one will be returned arbitrarily
    """
    needle = os.path.join(*parts)
    for match in glob.iglob(needle):
        return match
    return None


def glob_one_mult_roots(roots, *parts: str) -> Optional[str]:
    """
    We have data spread across /data/ and /data2/.
    Hopefully this is temporary, but for now, we use this function.
    """
    for root in roots:
        result = glob_one(root, *parts)
        if result is not None:
            return result
    return None


def get_agg_dir(agg_id: int) -> Optional[str]:
    """
    Locates an aggregation directory on our filesystem
    Current implementation relies on path naming convention
    """
    return glob_one_mult_roots(agg_roots, "LN*", f"aggregation-{agg_id}")


def get_seq_dir(flowcell_label) -> Optional[str]:
    """
    Locates an sequencing directory on our filesystem
    """
    return glob_one_mult_roots(seq_roots, f"*{flowcell_label}")


def get_flowcell_dir(label) -> Optional[str]:
    """
    Locates a flowcell directory on our filesystem
    """
    return glob_one_mult_roots(flowcell_roots, f"FC{label}*")


@app.command()
def flowcell(label: str):
    seq_dir = get_seq_dir(label)
    flowcell_dir = get_flowcell_dir(label)
    events = Events([])
    if seq_dir:
        events = events.merged(get_seq_events(seq_dir))
        #print("Sequencing dir:", seq_dir, get_elapsed_time(seq_dir))
        # events = get_seq_events(seq_dir)

        #print("Events:", events.to_tsv())
    if flowcell_dir:
        events = events.merged(get_flowcell_events(flowcell_dir))
        #print("Flowcell dir:  ", flowcell_dir, get_elapsed_time(flowcell_dir))
        #events = get_flowcell_events(flowcell_dir)
    #print(f"FC{label}")
    print(events.to_tsv())


@app.command()
def alignment(aln_id: int):
    pass

@app.command()
def trace(filename: str):
    events = Events(events=parse_trace_file(filename))
    print(events.to_tsv())

@app.command()
def aggregation(agg_id: int):
    aggdir = get_agg_dir(agg_id)
    if aggdir is None:
        log.error("Could not find agg directory")
        return
    trace = get_trace_contents(aggdir)
    events = Events(parse_trace(trace))
    print(events.to_tsv())
    #print("\n".join(str(e) for e in events))

    #print("Total real time", get_elapsed_time(aggdir))


def find_traces(rootpath: str):
    """ A generator that yields all traces in a directory """
    exclude_dirs = set("work")
    for root, dirs, files in os.walk(rootpath, topdown=True):
        dirs[:] = [d for d in dirs if d not in exclude_dirs]
        for x in files:
            if x == "trace.txt":
                yield os.path.join(root, x)
                # Stop descending in this directory
                dirs[:] = []
                break


@app.command()
def directory(root_dir: str):
    """ Crawls a directory reporting all traces found """

    first_trace = True
    for trace in find_traces(root_dir):
        events = parse_trace_file(trace)
        print(Events(events).to_tsv(print_header=first_trace).strip())
        first_trace = False


@app.command()
def parse_time(time: str):
    print(parse_trace_duration(time))


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell
# after importing without automatically running it
if __name__ == "__main__":
    app()
