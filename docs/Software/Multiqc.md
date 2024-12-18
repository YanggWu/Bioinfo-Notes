# MultiQC

**MultiQC** 是一个用于汇总生物信息学分析结果的工具，可以自动从指定目录中查找常见的生物信息学分析日志文件并生成一个整合报告，适用于对多个样本的质量控制或结果汇总。

<span class="ct">主要功能</span>

- **汇总分析结果**：支持多种生物信息学工具的输出文件。
- **生成 HTML 报告**：包含交互式图表和样本的比较。
- **数据导出**：支持多种格式的报告或数据导出。
- **模块化运行**：可选择性地启用或禁用某些工具的分析模块。

## 运行

直接指定要分析的文件路径
在指定的目录中查找所有支持的输出文件，并生成一个默认报告。

```sh
# 在当前目录中查找所有支持的输出文件，并生成一个默认报告。
multiqc .	# 默认html报告为： multiqc_report.html 报告数据目录：multiqc_data

# 生成自定义输出报告
multiqc . -o reports -n my_report.html
    
# 仅生成包含 fastqc 模块的报告。
multiqc . -m fastqc
```

`-n, --filename` 报告和文件前缀

`-o, --outdir` 指定输出目录

`-m, -module` 仅包含特定工具的模块

## 常用参数解析

=== "主选项"

    - `--force` / `-f`：覆盖已存在的报告文件。
    - `--config` / `-c`：指定要加载的配置文件。
    - `--cl-config`：直接在命令行中指定 YAML 格式的配置。
    - `--filename` / `-n`：自定义报告文件名，例如 `-n my_report.html`。
    - `--outdir` / `-o`：设置报告输出目录。
    - `--ignore` / `-x`：忽略某些分析文件（支持通配符）。
    - `--ignore-samples`：忽略指定的样本名称。
    - `--ignore-symlinks`：忽略符号链接的目录或文件。
    - `--file-list` / `-l`：提供包含文件路径的列表文件，MultiQC 将根据此列表扫描文件。

=== "模块选择"

    - `--module` / `-m`：指定运行的模块，可以多次使用，例如 `-m fastqc -m samtools`。
    - `--exclude` / `-e`：排除指定模块。

=== "样本命名"

    - `--dirs` / `-d`：在样本名称中添加目录路径。
    - `--dirs-depth` / `-dd`：控制路径深度，负值表示从路径的开头取值。
    - `--fullnames` / `-s`：保留完整的样本文件名作为样本名称。
    - `--fn_as_s_name`：使用日志文件名作为样本名称。
    - `--replace-names`：提供 TSV 文件，用于在生成报告时替换样本名称。

=== "报告定制"

    - `--title` / `-i`：设置报告标题。
    - `--comment` / `-b`：在报告顶部添加自定义注释。
    - `--template` / `-t`：指定报告模板（如 `default`, `simple`）。
    - `--sample-names`：提供 TSV 文件，在报告中更改样本按钮的显示名称。
    - `--sample-filters`：提供 TSV 文件，为报告添加样本显示/隐藏过滤器。
    - `--custom-css-file`：为报告添加自定义 CSS 文件。

=== "输出文件"

    - `--flat` / `-fp`：仅生成静态图像。
    - `--interactive` / `-ip`：仅生成交互式图表。
    - `--export` / `-p`：将图表另存为静态图像。
    - `--data-dir` / `--no-data-dir`：强制创建解析数据目录。
    - `--data-format` / `-k`：输出解析数据为其他格式（如 `tsv`、`csv`、`json` 或 `yaml`）。
    - `--zip-data-dir` / `-z`：压缩数据目录。
    - `--no-report`：不生成报告，仅导出数据和图表。
    - `--pdf`：生成 PDF 格式报告（需要安装 Pandoc）。

=== "运行行为"

    - `--verbose` / `-v`：增加日志输出的详细信息。
    - `--quiet` / `-q`：仅显示警告信息。
    - `--strict`：启用严格模式，不忽略任何异常。
    - `--development` / `--dev`：开发模式，导出未压缩的图表数据。
    - `--require-logs`：强制要求所有显式请求的模块都有日志文件。
    - `--profile-runtime`：在报告中添加运行时间分析。
    - `--profile-memory`：分析模块内存使用情况。
    - `--no-megaqc-upload`：禁止将报告上传到 MegaQC。

=== "清理选项"

    - `--clean-up` / `--no-clean-up`：在完成后移除临时目录和日志文件。
    - `--no-version-check`：禁用版本检查。