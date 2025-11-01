**Telescope+Encode数据**
- [论文和数据来源](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1006453)
- [TelescopeEncode Github](https://github.com/mlbendall/telescope?utm_source=chatgpt.com)
- 流程：
  - `parallel-fastq-dum`下载ENCODE数据
  - `flexbar`去接头
  - `bowtie2`+`samtools`比对+拼接
  - 直接从[Telescope官方仓库](https://github.com/mlbendall/telescope_annotation_db/tree/master/builds)获取配套的gtf注释
  - `telescope`计数

**featureCounts+TEtranscripts+Encode数据**
- 仍使用TelescopeEncode的数据
- [featureCounts介绍和使用](https://scienceparkstudygroup.github.io/ibed-bioinformatics-page/source/core_tools/featurecounts.html)
- [TEtranscripts Github](https://github.com/mhammell-laboratory/TEtranscripts)
- 流程：
  - 序列下载比对步骤同上
  - 下载RepeatMasker注释，并自己提取其中hERV的部分
  - `featureCounts`和`TEcount`两种方法计数
  - `TEtranscripts`对两组进行差异表达分析（每组可以有多个样本）

**TE_Transcript_Assembly**
- [pipeline代码](https://epigenome.wustl.edu/TE_Transcript_Assembly/tool.html)
- [相关论文](https://www.cell.com/cell/fulltext/S0092-8674(21)01104-1)
