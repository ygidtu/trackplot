# sashimi

## Theme
global

- items:
  - 字体大小
  - 字体
  - color map 
  - hspace
  - wspace
  - etc.
- functions:
  - blank(axis: mpl.axes.Axes)
  - withYAxis(axis: mpl.axes.Axes, **kwargs)
  - withXAxis(axis: mpl.axes.Axes, **kwargs)
  - withBothAxis(axis: mpl.axes.Axes, **kwargs)
  

## Genomic

## Input objects
- items:
  - path
  - alias
  - color
  - type: ["bam", "bw", "gtf", etc...]
- functions:
  - Check() # 检查文件类型和合法性 
  - Prepare()  # 检查index等
  - Load() # 根据不同文件要求读取范围结果
  
### GenomicRegion

所有涉及到基因组坐标的基础类

- items:
  - chromosome
  - start
  - end
  - strand
  - name
  - id
- property
  - __str__ -> f"{chrom}:{start}-{end}:{strand}"
  - __lt__
  - __gt__

### GenomicAnnotations(GenomicRegion)

- items
  - color
  - alias
  - type: Enum[Stoke, Focus, Line]

### region(GenomicRegion)

- items:
  - chromosome
  - start
  - end
  - strand
  - sequence: Optional[str]
  - indicator_lines: List[int]
  - strokes: List[]
  - focus: List[]
- property:
  - Label: -> 原sashimi图下方的注释
  - __dict__ -> graph coords
  - plot_strokes(axis: mpl.axes.Axes)
  - plot_focus(axis: mpl.axes.Axes)
  - plot_indicator(axis: mpl.axes.Axes)

## Framework

用于组织架构图像整体，主要用于记录字体大小、字体、dpi输出格式等参数以及作图

不绑定图片大小，matplotlib能够自适应调整输出图像的大小

- items:
  - plots: Dict{PlotInterface: IntRange}  # 图形以及index
  - theme: Theme  统一调整风格
- functions:
  - save(output: str, dpi: int=300)
  - add_track(type: str, **kwargs)  # 添加图形的，根据添加的类别，动态调gridspec的大小，分配index
    - 生成不同的obj: PlotInterface
  - plot(chrom: str, start: int, end: int)  # 统一调用各种参数和作图接口
  
  
## PlotInterface
- items:
  - alias
  - color
  - data
  - number # 要占用几个grid
- property
  - size <- 用于Framework动态调整图像大小
- functions:
  - Load(path: str, **kwargs): Load data
  - Plot(ax: mpl.Axes, **kwargs) -> idx # stepped index 下一幅图的位置，统统使用axes来作图，偏于控制图像的整体布局

### Sashimi(PlotInterface)

- items:
  - data: List[objects]
  - size: 1 if not show_size else 2
  - show_side: bool
  - depth
  - plus, minus
  - junction: {Region: number}

- functions:
  - Load(path: str, side_plot: bool, strand = ["+", "-", None])   # 加载
  - Plot(ax: mpl.Axes, **kwargs) -> idx # stepped index 下一幅图的位置
    - draw_density(ax_var):
    - draw_density_by_strand(ax_var, draw_hist: bool)

### Lines(PlotInterface)

- items:
  - data: List[objects]
  - size: 1
  - do_scale
- functions:
  - Load(path: str, region: Region)  # 从bam、bigwig中加载信号，做成线图
  - Plot(ax: mpl.Axes, **kwargs) -> idx # stepped index 下一幅图的位置
     - draw_line
     - draw_lines_at_same_plot

### Heatmap(PlotInterface)

- items:
  - data: List[objects]
  - clustering
    - method
    - metric
  - size: 1
  - do_scale
- functions:
  - Load(path: str, type = ["bigwig", "list", "depth"], region: Region) # 从一个bigwig，samtools depth等文件中加载信号矩阵 
  - Plot(ax: mpl.Axes, **kwargs) -> idx # stepped index 下一幅图的位置，数据转化为matrix并处理

### Transcripts:

- items:
  - PATH
  - alias
  - data: {Region: List[region]} # 保存 transcript: List[exons]
  - size: 1
  - filter_empty
- functions:
  - Load(path: str, region: Region) # 从gtf，gff加载转录本信息
  - Plot(ax: mpl.Axes, **kwargs) -> idx # stepped index 下一幅图的位置

