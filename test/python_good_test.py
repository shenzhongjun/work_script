#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
个人理解-OOP的作用
提高代码可读性，降低维护难度
使代码模块化，可重复使用，通过继承等特性减少代码重复
使代码精简
"""


class Fruit(object):

    def __init__(self, insale):
        self.name = 'fruit'
        self.price = '5$'
        self.weight = None
        if insale: print(f'to buy {self.name} by {self.price} is very good! come on!')

    def heavy(self):
        return f'good, I am {int(self.weight) + 1000}g now!'

    @classmethod  # 第一个参数必须是cls，可以调用本身的属性、方法
    def red(cls):
        return f'{cls.__name__} is red. and {cls.big}'

    @staticmethod
    def blue():  # 可以接受参数但不能调用本身的属性、方法
        return 'I am blue, don`t eat!'

    x = 'I am big, have a bear please.'

    @property  # 调用类的属性需要先实例化，如果不写@big.setter，则此属性是只读的。返回值依赖外部变量，才能被调用。
    def big(self):
        return self.x

    @big.setter
    def big(self, value):
        self.x = value


class Apple(Fruit):

    def __init__(self, insale):
        super().__init__(insale)
        print(self.name)
        self.name = 'Apple'
        self.weight = '200g'
        print('new name:', self.name, self.weight, sep='\n')


if __name__ == '__main__':
    fruit = Fruit(insale=True)
    apple = Apple(insale=True)
    juice = Fruit(insale=False)
    juice.name = 'juice'
    juice.weight = '100'
    juice.big = 'juice is not big but needle.'
    print(juice)
    fruit.x = ''

    with open('dd') as x, \
            open('zz', 'w') as kk:
        pass
