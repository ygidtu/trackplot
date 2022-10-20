from sashimi.plot_func import *

from matplotlib import pyplot as plt


def test_density():
    from sashimi.file.Bam import Bam
    fig, ax = plt.subplots()
    region = GenomicLoci("chr1", 1270656, 1284730, "+")
    bam = Bam.create("../example/bams/1.bam")
    bam.load(region, log_trans="2")
    plot_density(ax, bam)
    plt.savefig("plot_density_bam.png")


def test_bw():
    from sashimi.file.Bigwig import Bigwig

    bw = Bigwig.create("../example/bws/1.bw", title="test")
    bw.load(GenomicLoci("chr1", 1270656, 1284730, "+"))
    fig, ax = plt.subplots()
    plot_density(ax, bw)
    plt.savefig("plot_density_bw.png")


def test_ref():
    region = GenomicLoci("chr1", 1270656, 1284730, "+")
    fig, ax = plt.subplots()
    ref = Reference.create(
        "../example/example.sorted.gtf.gz",
        primary=True,
        add_domain=False)
    ref.load(region)

    ref.add_interval(
        interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",
        label="PolyASite"
    )

    ref.add_interval(
        interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.simple.bed.gz",
        label="PolyASite_simple"
    )

    plot_reference(ax, ref,
                   show_gene=True,
                   show_id=True,
                   plot_domain=True,
                   show_exon_id=True
                   )

    plt.savefig("plot_reference.pdf")


def test_depth():
    from sashimi.file.Depth import Depth
    region = GenomicLoci("chr1", 1270656, 1284730, "+")
    depth = Depth.create("../example/depth.bgz")
    depth.load(region)
    fig, ax = plt.subplots(nrows=len(depth))
    idx = 0
    for x, y in depth.items():
        plot_density(ax[idx], region=region, data=y, y_label=x,
                     theme="ticks_blank" if idx < len(depth) - 1 else "ticks")
        idx += 1

    plt.savefig("plot_density_depth.png")


def test_heatmap_and_line():
    from matplotlib import gridspec
    from sashimi.file.Bam import Bam

    region = GenomicLoci("chr1", 1270656, 1284730, "+")
    graph_coords = init_graph_coords(region)
    data = {}
    attrs = {}
    colors = ["red", "blue", "yellow", "green", "grey"]
    for i in range(1, 5):
        bam = Bam.create(f"../example/bams/{i}.bam")
        bam.load(region)
        data[str(i)] = bam.data
        attrs[str(i)] = {"color": colors[i]}

    gs = gridspec.GridSpec(1, 2, width_ratios=(.99, .01), wspace=0.01, hspace=.15)
    ax_var = plt.subplot(gs[:, 0]),
    cbar_ax = plt.subplot(gs[:, 1])
    plot_heatmap(ax_var[0], cbar_ax, data, do_scale=True, clustering=True)
    plt.savefig("plot_heatmap.png")

    fig, ax = plt.subplots()
    plot_line(ax, data, line_attrs=attrs, graph_coords=graph_coords)
    plt.savefig("plot_line.png")


def test_graph_coord():
    from matplotlib import gridspec
    region = GenomicLoci(chromosome="1", start=100, end=800, strand="+")
    exon_width = .5
    exons = [
        [150, 300], [320, 450],
        [460, 500], [470, 500],
        [470, 600], [700, 800],
    ]

    gs = gridspec.GridSpec(2, 1)
    ax = plt.subplot(gs[0, 0])
    for ind, exon in enumerate(exons):
        s, e, strand = exon[0], exon[1], "+"
        x = [s, e, e, s]
        x = [i - region.start for i in x]
        y = [
            ind - exon_width / 2, ind - exon_width / 2,
            ind + exon_width / 2, ind + exon_width / 2
        ]
        ax.fill(x, y, 'k', lw=.5, zorder=20)

    ax.set_xbound(0, len(region))
    graph_coords = init_graph_coords(region, exons, intron_scale=.5)

    ax = plt.subplot(gs[1, 0])
    for ind, exon in enumerate(exons):
        s, e, strand = exon[0] - region.start, exon[1] - region.start, "+"
        x = [graph_coords[s], graph_coords[e], graph_coords[e], graph_coords[s]]
        y = [
            ind - exon_width / 2, ind - exon_width / 2,
            ind + exon_width / 2, ind + exon_width / 2
        ]
        ax.fill(x, y, 'k', lw=.5, zorder=20)
    ax.set_xbound(0, len(region))
    plt.savefig("plot_graph_coords.png")


def test_set_x_ticks():
    region = GenomicLoci(chromosome="1", start=100, end=120, strand="+")
    fig, ax = plt.subplots()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False)
    set_x_ticks(ax, region, font_size=10, sequence={110: "A", 106: "C"})
    plt.savefig("plot_x_ticks.png")


def test_igv_plot():
    from sashimi.file.ReadSegments import ReadSegment
    fig, ax = plt.subplots()
    region = GenomicLoci("chr1", 13362, 29900, "+")

    rs = ReadSegment.create("../example/bams/WASH7P.bam",
                            features={"m6a": "ma", "real_strand": "rs", "polya": "pa"})
    rs.load(region)
    plot_igv_like(ax, {'full': rs}, y_label="fl")
    plt.savefig("test_igv_plot.pdf")


def test_igv_plot2():
    from sashimi.file.ReadSegments import ReadSegment
    fig, ax = plt.subplots()
    region = GenomicLoci("chr1", 1270656, 1284730, "+")

    rs = ReadSegment.create("../example/bams/0.bam",
                            features={"m6a": "ma", "real_strand": "rs", "polya": "pa"})
    rs.load(region)
    plot_igv_like(ax, {'full': rs}, y_label="fl")
    plt.savefig("test_igv_plot.2.pdf")


def test_igv_plot3():
    from sashimi.file.ReadSegments import ReadSegment
    fig, ax = plt.subplots()
    # 1: 10024601 - 10038168
    region = GenomicLoci("1", 10024601, 10038168, "+")
    rs = ReadSegment.create("../example/test.bed.gz")
    rs.load(region)
    plot_igv_like(ax, {'full': rs}, y_label="fl")
    plt.savefig("test_igv_plot.3.pdf")


def test_hic_plot():
    from sashimi.file.HiCMatrixTrack import HiCTrack
    fig, ax = plt.subplots()
    region = GenomicLoci("X", 2500000, 2600000, "*")
    hic = HiCTrack.create(path="../example/Li_et_al_2015.h5", label="Li",
                          depth=30000,
                          trans="log2"
                          )
    hic.load(
        region=region
    )
    plot_hic(ax, {hic.label: hic}, y_label=hic.label)
    plt.savefig("test_hic.pdf")


def test_motif():
    fig, ax = plt.subplots()

    motif = Motif.create("../example/motif.bed.gz", GenomicLoci("chr1", 100, 105, "*"))
    motif.load(GenomicLoci("chr1", 100, 105, "*"))
    print(motif.data)

    plot_motif(ax, motif, width=0.8)
    plt.savefig("test_motif.pdf")


# test_igv_plot()
# test_igv_plot2()
# test_ref()
test_motif()
