## 目标
- 当 `fftN > 1024` 时自动开启 `usePowerSpectrum`，降低高 FFT 点数下的计算开销，同时保持对外行为稳定。

## 核心改动
1. 在 `src/funkin/vis/dsp/SpectralAnalyzer.hx` 增加开关
- 新增 `public var usePowerSpectrum(default, set):Bool = false`。
- 在 `set_fftN(value:Int)` 中自动切换：`usePowerSpectrum = (pow2 > 1024)`；当 `pow2 <= 1024` 时关闭。

2. 在 `git/src/grig/audio/FFT.hx` 增加功率谱计算
- 新增 `calcFreqPower(data:Array<Float>):Array<Float>`：使用 `real^2 + imag^2` 输出功率谱，复用内部工作缓冲；保留 `calcFreq()` 作为幅度谱（兼容旧逻辑）。

3. 在 `SpectralAnalyzer.getLevels()` 选择路径
- 当 `usePowerSpectrum` 为真：调用 `fft.calcFreqPower(signal)` 得到功率谱。
- 传入可视化聚合时，用新参数标记功率谱路径（见下一步）。

4. 在 `git/src/grig/audio/FFTVisualization.hx` 支持功率谱
- 扩展 `makeLogGraph(freq:Array<Float>, bands:Int, dbRange:Int, intRange:Int, usePower:Bool = false)`。
- 传递到 `computeFreqBand(...)`，将 dB 映射改为 `k * log10(n)`，其中 `k = usePower ? 10 : 20`。
- 保留已实现的动态基数与边界保护（`base = freq.length`）。

## 自动启用与兼容
- 默认行为保持不变（幅度谱）。仅当 `fftN` 被设置为 >1024（或构造为此值）时自动切换到功率谱路径。
- 若需要强制某一模式，可通过 `usePowerSpectrum` 显式设置覆盖自动切换。

## 验证与调参
- 在 `fftN=2048`、`barCount=128`、`updateIntervalMs=33` 场景下对比 CPU 时间，预期有可测降幅。
- 验证条带高度与趋势与旧版一致性（允许常量因子差异），必要时微调 `dbRange` 或条带加权。

## 交付
- 修改点集中在上述四处，接口保持向后兼容；新增一个可选参数与一个布尔开关。请确认后我将开始实施并提交对应代码修改。