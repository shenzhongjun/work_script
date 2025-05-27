conda activate pytorch

# 将视频转为音频
ffmpeg -i S40323-21304266.mp4 20240323课.mp3

# 音频降噪
ffmpeg -i input/22.1.29.mp3 -af "volume=4,anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" output/22.1.29_new.aac
ffmpeg -i input/20221026.m4a -af "anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" output/20221026.m4a

# 批量降噪 
for i in `ls input/*`;do name=`echo $i|cut -f2 -d/`;echo $i; ffmpeg -i $i -threads 4 -af "volume=4,anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" output/$name; done

# 音频合并左右声道，aac转为mp3效果比保持aac好 课
ffmpeg -i 20231026.aac -ac 1 20231026.mp3

# 批量合并左右声道
for i in `ls input/*`;do name=`echo $i|cut -f2 -d/`;echo $i; ffmpeg -i $i -ac 1 output/$name; done

# 音频1和音频2拼接
ffmpeg -i "concat:input/1.mp3|input/2.mp3" -acodec copy output/output.mp3

# ----------测试用命令-----------
# 将视频转为音频。实测转为aac拖不动，没法听！
ffmpeg -i 22.1.29.mp4 22.1.29.aac
# 将视频转为音频并降噪。实测转完后没法听！
ffmpeg -i 22.1.29.mp4 -threads 4 -af "volume=4,anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" test.aac
# 音频降噪
# -threads 4 当前版本6不支持多线程，7可能会支持
# -af：表示应用音频滤镜
# "volume=4": 音量调为4倍
# "anlmdn"：是自适应降噪滤镜，用于减小背景噪音，并尽量保留人声。
# "highpass=f=100,lowpass=f=3000"：这是一个带通滤波器，用于滤除频率低于100Hz和高于3000Hz的噪声。
# "equalizer=f=1000:t=h:w=150:g=5"：这是一个均衡器滤镜，用于增强1000Hz频率附近的声音，并减小其他频率的声音。
# "acompressor=threshold=-20dB:ratio=4:attack=50:release=100"：这是一个动态范围压缩器滤镜，用于平衡音频的动态范围。其中threshold 设置了压缩的阈值，ratio 设置了压缩的比率，attack 和 release 分别设置了压缩器的攻击时间和释放时间。
ffmpeg -i 22.1.29.mp3 -threads 4 -af "volume=4,anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" 22.1.29_new.aac
# 对于本来音量就够的音频无需提高音量，反而会增加一些噪音！比特率也不要随便增加，会导致音频拖不动卡住！
ffmpeg -i input/20221026.m4a -af "anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" output/20221026.m4a
# 批量降噪
for i in `ls input/*`;do name=`echo $i|cut -f2 -d/`;echo $i; ffmpeg -i $i -threads 4 -af "volume=4,anlmdn,highpass=f=100,lowpass=f=3000,equalizer=f=1000:t=h:w=150:g=5,acompressor=threshold=-20dB:ratio=4:attack=50:release=100" output/$name; done
# 音频合并左右声道-用于小米内录，以下实测皆可用，但是音轨和时长会对不上，不知道怎么解决！先用第一行的！
ffmpeg -i input.mp3 -af "pan=stereo|c0=c0+c1|c1=c0+c1" output.mp3
ffmpeg -i "20240304三神合一.mp3" -af "pan=stereo|c0=c0+c1|c1=c0+c1" -c:a aac "20240304三神合一2.m4a"
ffmpeg -i input_file.mp3 -af "pan=stereo|c0=c0|c1=c0" -c:a libmp3lame output_file.mp3
ffmpeg -i 20240325-问阳神.aac -ac 1 20240325问阳神2.aac
# 批量合并左右声道，第一条命令加-b:a 192k速度慢且音频拖不动，拖到后面就没声音了！因此使用最后一条命令
for i in `ls input/*`;do name=`echo $i|cut -f2 -d/`;echo $i; ffmpeg -i $i -af "pan=stereo|c0=c0+c1|c1=c0+c1" -b:a 192k output/$name; done
ffmpeg -i input/20210718课.aac -ac 1 test.aac	# 实测总长度和拖动没问题，听起来也没问题，比特率默认降到mono, fltp, 69
ffmpeg -i input/20210718课.aac -ac 1 -b:a 192k test2.aac	# 实测最高96，拖动变慢，听起来并没有更清晰
for i in `ls input/*`;do name=`echo $i|cut -f2 -d/`;echo $i; ffmpeg -i $i -ac 1 output/$name; done