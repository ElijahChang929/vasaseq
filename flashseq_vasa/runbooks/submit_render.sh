#!/usr/bin/env bash
# submit_render.sh -- 把 render_runbook_direct.sh 作为 SLURM 作业提交到 nint。
#
# 为什么要有这一层（而不是直接敲 sbatch --wrap）：
#
#   1. 绝不能在登录节点渲染。登录节点的 memory cgroup 会 *静默* 杀死进程，
#      你拿到的是一个被截断的 HTML 和一个愉快的 exit 0。
#   2. `srun --pty` 会随着 SSH 断线一起死。VS Code 的终端断得很随意，
#      所以默认走 sbatch，而不是交互式。
#   3. 资源参数不是拍脑袋来的，见下面 SBATCH 块里的注释。
#
# 用法:
#   ./submit_render.sh                      # 渲染 full_project.qmd
#   ./submit_render.sh phase_rrna_readlength.qmd
#
# 提交后:
#   squeue -u $USER -n render-runbook
#   tail -f render-<jobid>.log              # Ctrl-C 只停止 tail，作业继续跑

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
QMD="${1:-full_project.qmd}"

# 相对路径按本脚本所在目录解析，这样从任何地方调用都对
[[ "$QMD" = /* ]] || QMD="$HERE/$QMD"
[ -f "$QMD" ] || { echo "找不到文件: $QMD" >&2; exit 2; }

[ -x "$HERE/render_runbook_direct.sh" ] || {
    echo "render_runbook_direct.sh 不存在或不可执行: $HERE" >&2; exit 2; }

# 分区必须是 ncpu，不能是 nint。实测 2026-07-30：
#   sbatch --partition=nint  ->  slurm_job_submit: refusing non-interactive job
#                                in interactive partition.
# nint 只接受交互作业（srun --pty）。批量渲染走 ncpu（默认分区，7 天上限）。
# nint 仍然是逐 chunk 调试的正确去处 —— 见 start_jupyter_nint.sh。
#
#   --cpus-per-task=4  quarto/pandoc 基本单线程，4 核是给 chunk 里可能的
#                      并行留的余量
#   --mem=16G          实测 MaxRSS = 716812K ≈ 0.7 GB (job 51055258, cn 节点)。
#                      16G 是 20 倍余量，够任何 chunk 的临时膨胀。
#   --time=01:00:00    实测 full_project.qmd = 00:04:12 (job 51055258)，
#                      其中约 40s 是首次 shell 冷启动预热。1 小时足够。
#                      活跃 chunk 只读日志和表；重活（demux/mapping）在文档里
#                      是不执行的 fenced block。
jid=$(sbatch --parsable \
    --partition=ncpu \
    --cpus-per-task=4 \
    --mem=16G \
    --time=01:00:00 \
    --job-name=render-runbook \
    --output="$HERE/render-%j.log" \
    --chdir="$HERE" \
    --wrap="$HERE/render_runbook_direct.sh '$QMD'")

echo "已提交作业 : $jid"
echo "渲染文件   : $QMD"
echo "日志       : $HERE/render-${jid}.log"
echo
echo "查看进度   : tail -f $HERE/render-${jid}.log"
echo "查看队列   : squeue -u $USER -n render-runbook"
echo
echo "作业结束后必须读日志末尾的 VERIFY 块 —— quarto 会把出错 chunk 的报错"
echo "内嵌进 HTML 然后照样 exit 0，退出码本身不能证明任何事情。"
