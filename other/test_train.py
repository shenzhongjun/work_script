# 导入库
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split

# 生成模拟数据，有1000个样本，20个特征（产出X：1000行20列的数据，即数据集），分为2个类别（产出y，每个样本的类别，即标签集）。
# n_informative默认为2，即每个类别有2个特征是有意义的（具有高区分度discriminatory power）。随机数设为42
X, y = make_classification(n_samples=1000, n_features=20, n_informative=2, n_classes=2, random_state=42)
# 把X和y转为torch张量
X = torch.tensor(X, dtype=torch.float32)
y = torch.tensor(y, dtype=torch.long)

# 划分数据集，训练集占0.8，测试集占0.2。大模型可以数据集占比高一些，但小数据集为了保证测试结果准确，最好占20%或更高
# 通常认为，样本数量超过10万个的数据集可以被称为大数据集，而小于10万个的数据集则被视为小数据集。
# 特征数量超过几百或几千个的数据集可以被视为大数据集，而特征数量较少的数据集则可以被视为小数据集。
# 200个样本，80万个特征的数据可以被视为高维稀疏数据，而不是大数据集。虽然特征数量非常大，但是样本数量非常少，可能会导致过拟合等问题。
# 在这种情况下，需要进行特殊的数据处理和特征选择，以提高模型的泛化能力和性能。同时，由于数据集的维度非常高，
# 可能需要使用一些高维数据处理和降维技术，如主成分分析（PCA）和 t-SNE 等，来降低数据集的维度，以便进行可视化和分析。
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# 定义模型
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(500, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, 2)
        self.relu = nn.ReLU()
        self.softmax = nn.Softmax(dim=1)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.relu(self.fc2(x))
        x = self.softmax(self.fc3(x))
        return x


# 实例化模型和优化器
net = Net()
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(net.parameters(), lr=0.001)

# 训练模型
for epoch in range(20):
    running_loss = 0.0
    for i in range(len(X_train)):
        optimizer.zero_grad()
        outputs = net(X_train[i])
        loss = criterion(outputs.unsqueeze(0), y_train[i].unsqueeze(0))
        loss.backward()
        optimizer.step()
        running_loss += loss.item()
    print('Epoch %d loss: %.3f' % (epoch + 1, running_loss / len(X_train)))

# 在测试集上进行预测
outputs = net(X_test)
_, y_pred = torch.max(outputs, 1)

# 如果这样一个高维稀疏数据使用t-SNE降维后使用pytorch进行深度学习，怎么操作，请给出具体代码，谢谢
# t-SNE 是一种非线性降维技术，可以用于高维数据可视化和特征提取。在使用 t-SNE 对高维稀疏数据进行降维之后，可以使用 PyTorch 构建深度学习模型，例如全连接神经网络或卷积神经网络等。下面是一个简单的示例代码，演示了如何使用 t-SNE 进行降维和 PyTorch 进行分类。

# 首先，我们可以使用 Scikit-learn 的 t-SNE 函数对数据进行降维：
import numpy as np
from sklearn.manifold import TSNE

# 读取高维稀疏数据
X = np.load('data.npy')

# 使用 t-SNE 进行降维
X_embedded = TSNE(n_components=2).fit_transform(X)

# 然后，我们可以使用 PyTorch 构建一个简单的全连接神经网络，并对数据进行分类：
import torch
import torch.nn as nn
import torch.optim as optim


# 定义全连接神经网络
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(2, 10)
        self.fc2 = nn.Linear(10, 5)
        self.fc3 = nn.Linear(5, 2)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x


# 将数据转换为 PyTorch 张量
X_embedded = torch.Tensor(X_embedded)
y = torch.Tensor(np.random.randint(2, size=X_embedded.shape[0]))

# 定义损失函数和优化器
criterion = nn.CrossEntropyLoss()
optimizer = optim.SGD(net.parameters(), lr=0.01)

# 训练模型
net = Net()
for epoch in range(100):
    optimizer.zero_grad()
    output = net(X_embedded)
    loss = criterion(output, y.long())
    loss.backward()
    optimizer.step()

# 使用模型进行预测
predicted = torch.argmax(net(X_embedded), dim=1)


# 请问上面这个数据集样本量是不是太少了，神经网络分类能取得较好的效果吗？如果换成传统机器学习分类方法会怎么样？如果要增加样本获得更好的性能，样本数量加到多少合适呢？
# 是的，上面这个数据集的样本量较少，如果直接使用神经网络可能会过拟合，导致在测试集上的表现不好。相比之下，传统机器学习分类方法可能会更适合这种数据集，因为它们通常需要较少的数据量就可以取得良好的效果。
# 如果想要在这个数据集上使用神经网络获得更好的性能，需要增加样本数量。具体来说，增加的样本数量应该是有一定规律的，即需要考虑样本数量与特征数量的比例，同时需要考虑模型复杂度的影响。在增加样本数量时，可以采用数据增强等方法来扩充数据集。具体来说，可以使用图像处理库对图像进行随机旋转、裁剪、缩放等变换，以生成更多的训练样本。此外，还可以使用合成数据的方法，例如在合成图像上添加噪声、改变亮度、对比度等来增加训练数据的多样性。

# 请注意，我们是在生物信息学背景下讨论神经网络的应用，前面的例子里也是假设神经网络用于基因表达数据的分析，不存在使用图像处理库对图像进行变换的可能
# 非常抱歉，我的回答偏离了生物信息学的应用场景。在生物信息学领域，数据集的样本量通常会受限于实验设计、数据采集的成本、时间等因素。在这种情况下，为了避免过拟合，可以采用一些措施来减轻模型复杂度，例如增加正则化项，降低网络层数，增加dropout等。同时，可以尝试使用其他机器学习方法，如随机森林、支持向量机等，它们的泛化能力通常会更好一些。如果数据集较小，可以考虑使用交叉验证等方法来评估模型的性能。



# 学习率是指模型在每次更新参数时，沿着梯度方向前进的步长1 。学习率太大，可能导致模型跳过最优解或者震荡不收敛；学习率太小，可能导致模型收敛速度太慢或者陷入局部最优1 。
# 正则化参数是指模型在损失函数中加入的一些额外的项，用来惩罚模型的复杂度，防止过拟合2 。正则化参数太大，可能导致模型欠拟合或者损失信息；正则化参数太小，可能导致模型过拟合或者过于敏感2 。
# 网络结构是指模型的基本架构，例如使用哪些类型的神经网络层，以及它们之间的连接方式2 。网络结构决定了模型的表达能力和适应性，不同的网络结构可能适合不同的任务和数据集2 。
# 批次大小是指模型在每次训练时，从数据集中抽取的样本数量2 。批次大小太大，可能导致模型训练速度慢或者内存不足；批次大小太小，可能导致模型训练不稳定或者收敛困难2 。









