# cython: profile=False

import click
from multiprocessing import cpu_count
from importlib.metadata import version
from dodi import input_stream_alignments
import logging
import datetime
import time


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
dodi_version = version("dodi")


@click.command()
@click.argument("sam", type=str, default='-', required=False)
@click.argument("output", required=False, type=click.Path())
@click.option("--paired", help="Paired end reads (or single)", default=defaults["paired"],
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-I', '--template-size', help="Insert size and stdev as 'FLOAT,FLOAT'",
              default=defaults["template_size"], type=str, show_default=False)
@click.option("-K", "--keep-mapq", help="Keep original mapq (applies to paired-end only)", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-S", "--secondary", help="Output secondary mappings", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--min-aln", help="Minimum alignment length", default=defaults["min_aln"], type=float, show_default=True)
@click.option("--max-overlap", help="Maximum overlap between successive alignments", default=defaults["max_overlap"],
              type=float, show_default=True)
@click.option("--ins-cost", help="Insertion cost", default=defaults["ins_cost"], type=float, show_default=True)
@click.option("--ol-cost", help="Overlapping alignment cost", default=defaults["ol_cost"], type=float,
              show_default=True)
@click.option("-c", "--inter-cost", help="Cost of inter-chromosomal jump", default=defaults["inter_cost"], type=float,
              show_default=True)
@click.option("-u", '--pairing-cost', help="Pairing cost", default=defaults["u"], type=float, show_default=True)
@click.option("-m", "--match-score", help="Matched base score", default=defaults["match_score"],
              type=float, show_default=True)
@click.option('--include', help=".bed file, elevate alignment scores in these regions. Determined by '--bias'",
              default=None, type=click.Path(exists=True))
@click.option("--bias", help="""Multiply match score by bias if alignment falls within regions .bed file
Unused if .bed not provided""", default=defaults["bias"], type=float, show_default=True)
@click.option("--logging", help="Set logging level", default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING']), show_default=False)
@click.version_option()
def dodi_choose(**kwargs):
    """Choose an optimal set of alignments from a collection of candidates.
    If reads are paired, alignments must be sorted by read-name with the bit flag
    designating read_1 vs read_2."""
    t0 = time.time()
    if kwargs['logging'] == 'INFO':
        rootLogger.setLevel(logging.INFO)
    elif kwargs['logging'] == 'DEBUG':
        rootLogger.setLevel(logging.DEBUG)
    elif kwargs['logging'] == 'WARNING':
        rootLogger.setLevel(logging.WARNING)
    logging.info(f"dodi v{dodi_version}")
    logging.info(str(kwargs))
    kwargs['u'] = kwargs['pairing_cost']
    kwargs['paired'] = kwargs['paired'] == 'True'
    kwargs['modify_mapq'] = not kwargs['keep_mapq']
    kwargs['version'] = dodi_version
    if kwargs['template_size'] != 'auto':
        insert_std = kwargs["template_size"].split(",")
        kwargs["insert_median"] = float(insert_std[0])
        kwargs["insert_stdev"] = float(insert_std[1])
    else:
        kwargs["insert_median"] = 300.
        kwargs["insert_stdev"] = 200.

    if not kwargs["include"]:
        kwargs["bias"] = 1.0
    else:
        logging.info("Elevating alignments in --include with --bias {}".format(kwargs["bias"]))

    if kwargs['paired']:
        batch_size = 10_000
    else:
        batch_size = 100

    count = input_stream_alignments.process_reads(kwargs, batch_size)
    logging.info("dodi processed {} query names in {} h:m:s".format(count, str(datetime.timedelta(
        seconds=int(time.time() - t0)))),
                 )


def test_command(ctx, **kwargs):
    """Run dysgu tests"""
    pass
    # tests_path = os.path.dirname(__file__) + "/tests"
    #
    # runner = CliRunner()
    # t = [tests_path + '/small.dysgu.sam']
    # click.echo(t)
    # result = runner.invoke(dodi_aligner, t)

