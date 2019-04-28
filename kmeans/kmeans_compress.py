"""
一张彩色图的大小是：（长*宽*3），因为每一个像素点都有 3 个颜色通道，相当于
对于每一个像素点，它的红黄蓝三种颜色都各有一个在 [0, 255] 区间内的值
"""

from skimage import io
from sklearn.cluster import KMeans
import numpy as np

image = io.imread('rabbit.jpg')
#io.imshow(image)
#io.show()

rows = image.shape[0]
cols = image.shape[1]

image = image.reshape(image.shape[0] * image.shape[1], 3)  # 3 代表三个颜色通道，是彩色图
# reshape 之后，变成 n 行 3 列的矩阵，每一行都是一个样本点，分别是各自红黄蓝三种颜色通道的值
kmeans = KMeans(n_clusters=128, n_init=10, max_iter=200)  # 分成 128 个簇
kmeans.fit(image)

clusters = np.asarray(kmeans.cluster_centers_, dtype=np.uint8)
labels = np.asarray(kmeans.labels_, dtype=np.uint8)  # labels_指的是每个样本点所属簇的 id ，取值从 0-128
labels = labels.reshape(rows, cols)
"""
到此为止，图片已经压缩完，原始图片大小为[长，宽，(红，黄，蓝)]，
压缩后图片大小为[长，宽，cluster id]，
所以由彩色图压缩成了灰度图
"""

print(clusters.shape)  # (128, 3)
io.imsave('compressed_rabbit.jpg', labels)

image = io.imread('compressed_rabbit.jpg')
io.imshow(image)
io.show()
