# cython: profile=False

import click
from multiprocessing import cpu_count
import pkg_resources
from . import input_stream_alignments
import logging
import sys
from sys import stderr


logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.WARNING)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)


defaults = {"paired": "True",
            "template_size": "auto",
            "replace_hardclips": "False",
            "max_insertion": 150,
            "min_aln": 10,

            "max_overlap": 5000,
            "ins_cost": 0.1,
            "ol_cost": 2,

            "zero_cost_boundary": 0,
            "max_gap_cost": 100,

            "inter_cost": 1,
            "u": 9,
            "match_score": 1,
            "bias": 1.15,
            "secondary": 'True'
            }

cpu_range = click.IntRange(min=1, max=cpu_count())
version = pkg_resources.require("dodi")[0].version


def show_params():
    args = " ".join(sys.argv[1:])
    logging.info(f"[dodi] Version: {version} {args}")


@click.command()
@click.argument("sam", type=str, default='-', required=False)
@click.argument("output", required=False, type=click.Path())
@click.option("--paired", help="Paired end reads (or single)", default=defaults["paired"],
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-I', '--template-size', help="Insert size and stdev as 'FLOAT,FLOAT'",
              default=defaults["template_size"], type=str, show_default=False)
# @click.option('-Q', '--modify-mapq', help="For paired-end reads, modify split-read MapQ ", default='True',
#               type=click.Choice(['True', 'False']), show_default=True, required=False)
# @click.option("-S", "--secondary", help="Output secondary mappings", type=click.Choice(["True", "False"]),
#               show_default=True, default=defaults["secondary"])
@click.option("-K", "--keep-mapq", help="Keep original mapq (applies to paired-end only)", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-S", "--secondary", help="Output secondary mappings", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-T", "--tags", help="Add sam tags to output", is_flag=True, flag_value=True, show_default=False, default=False)
# @click.option("--replace-hardclips",  help="Replace hard-clips with soft-clips when possible",
#               default=defaults["replace_hardclips"], type=click.Choice(["True", "False"]), show_default=True)
# @click.option("--tags",  help="Add sam tags to output",
#               default="False", type=click.Choice(["True", "False"]), show_default=True)
# @click.option("--fq1",  help="Fastq/fasta reads 1, used to add soft-clips to all hard-clipped read 1 alignments",
#               default=None, type=click.Path(), show_default=True)
# @click.option("--fq2",  help="Fastq/fasta reads 2, used to add soft-clips to all hard-clipped read 2 alignments",
#               default=None, type=click.Path(), show_default=True)
# @click.option("--max-insertion", help="Maximum insertion within read", default=defaults["max_insertion"], type=float,
#               show_default=True)
@click.option("--min-aln", help="Minimum alignment length", default=defaults["min_aln"], type=float, show_default=True)
@click.option("--max-overlap", help="Maximum overlap between successive alignments", default=defaults["max_overlap"],
              type=float, show_default=True)
@click.option("--ins-cost", help="Insertion cost", default=defaults["ins_cost"], type=float, show_default=True)
@click.option("--ol-cost", help="Overlapping alignment cost", default=defaults["ol_cost"], type=float,
              show_default=True)
# @click.option("-z", "--zero-cost-boundary", help="Gaps < z in length receive zero cost while gaps >= z have a linear cost", default=defaults["zero_cost_boundary"], type=float, show_default=True)
# @click.option("-g", "--max-gap-cost", help="Gaps >= g in length receive max-gap-cost", default=defaults["max_gap_cost"], type=float,
#               show_default=True)
@click.option("-c", "--inter-cost", help="Cost of inter-chromosomal jump", default=defaults["inter_cost"], type=float,
              show_default=True)
@click.option("-u", '--pairing-cost', help="Pairing cost", default=defaults["u"], type=float, show_default=True)
@click.option("-m", "--match-score", help="Matched base score", default=defaults["match_score"],
              type=float, show_default=True)
# @click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--include', help=".bed file, elevate alignment scores in these regions. Determined by '--bias'",
              default=None, type=click.Path(exists=True))
@click.option("--bias", help="""Multiply match score by bias if alignment falls within regions .bed file
Unused if .bed not provided""", default=defaults["bias"], type=float, show_default=True)
@click.option("--logging", help="Set logging level", default='WARNING',
              type=click.Choice(['INFO', 'WARNING']), show_default=False)
def dodi_aligner(**kwargs):
    """Choose an optimal set of alignments from a collection of candidates.
    If reads are paired, alignments must be sorted by read-name with the bit flag
    designating read_1 vs read_2."""

    if kwargs['logging'] == 'INFO':
        rootLogger.setLevel(logging.NOTSET)

    show_params()
    print('dodi parameters', kwargs, file=sys.stderr)
    kwargs['u'] = kwargs['pairing_cost']
    kwargs['paired'] = kwargs['paired'] == 'True'
    kwargs['modify_mapq'] = not kwargs['keep_mapq']
    kwargs['procs'] = 1

    input_stream_alignments.process_reads(kwargs)


def test_command(ctx, **kwargs):
    """Run dysgu tests"""
    pass
    # tests_path = os.path.dirname(__file__) + "/tests"
    #
    # runner = CliRunner()
    # t = [tests_path + '/small.dysgu.sam']
    # click.echo(t)
    # result = runner.invoke(dodi_aligner, t)

