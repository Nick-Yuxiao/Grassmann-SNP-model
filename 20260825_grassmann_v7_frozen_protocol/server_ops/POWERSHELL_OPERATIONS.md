# PowerShell 操作顺序

以下命令由用户在 Windows PowerShell 中执行。脚本默认 `Plan`，仅显示路径，
不会连接服务器。密码或私钥口令不要写入脚本参数、日志或聊天。

```powershell
$ops = 'C:\Users\Yuxiao Tan\OneDrive\桌面\Qinghua_bioinfo\Grassmann_model\deliverables\20260825_grassmann_v7_frozen_protocol\server_ops\GrassmannServerOps.ps1'

# 0. 本地预览，不连接服务器
& $ops -Action Plan -ReleaseId 'v7_20260825_p0'

# 1. 上传到新的、不可覆盖的 release 目录，并远端校验 SHA-256
& $ops -Action UploadV7 -ReleaseId 'v7_20260825_p0'

# 2. 建立 v6/v7 分离目录并保存第一次资源/任务快照
& $ops -Action Bootstrap -ReleaseId 'v7_20260825_p0'

# 3. 单独执行只读审计；这一步不会占 GPU
& $ops -Action Audit -ReleaseId 'v7_20260825_p0'

# 4. 建立隔离 cu128 环境并做最小 CUDA smoke；无空闲 GPU 时自动退出
& $ops -Action T00 -ReleaseId 'v7_20260825_p0'

# 5. RTX 5090 矩阵前向/反向测试；占用变化时自动退出
& $ops -Action GpuTest -ReleaseId 'v7_20260825_p0'

# 6. GPU 测试通过后才启动 T03 100-step profiler
& $ops -Action T03 -ReleaseId 'v7_20260825_p0' -T03OutputId 't03_profile_20260825'
```

如果当前机器禁止直接运行 `.ps1`，不要修改系统全局策略；可只对这一进程使用：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File $ops -Action Plan -ReleaseId 'v7_20260825_p0'
```

每一步是独立动作。不要一次粘贴整段执行；查看上一步退出码和服务器上生成的
结果后再继续。若 SSH 公钥仍未恢复，PowerShell 会在第一个远端动作停止，不会
改变服务器文件。

## 结果位置

- v6 结果：`/data1/home/tanyuxiao/Grassmann_model/v6/results/`
- v6 资源：`/data1/home/tanyuxiao/Grassmann_model/v6/resources/`
- v7 结果：`/data1/home/tanyuxiao/Grassmann_model/v7/results/`
- v7 资源：`/data1/home/tanyuxiao/Grassmann_model/v7/resources/`
- v7 release：`/data1/home/tanyuxiao/Grassmann_model/v7/code/releases/v7_20260825_p0/`

所有测试输出均使用新时间戳目录；脚本没有删除、终止或抢占命令。
