# kmeans用于图像压缩

## 环境
Windows7 + pycharm + anaconda3

## 文件清单
* kmeans_compress.py
* rabbit.jpg

## 代码说明
1. 用到的库有`skimage, sklearn.cluster, numpy`
2. 输入文件为图片，输出文件为压缩后的图片
3. `imread`函数返回的是`ndarray`，行和列分别表示图片的长和宽，即800*800；每个元素又是一个三维向量，分别表示该像素的RGB
4. `reshape`的意义是将图片拉直，转化为一维数组。由于是彩色图，每个像素点有三个颜色通道
5. KMeans函数：


