#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
测试内容：
1. @setter装饰器：GPT解释：动态更新？？？
2. @property装饰器：测试成功，可以实现动态更新！即每次调用类的示例的该属性，都会调用一次@property装饰的函数！

"""


class Cart:
    def __init__(self):
        self.items = []
        # self.total_price = 0

    def add_item(self, item):
        self.items.append(item)
        # self.calculate_total_price()

    def remove_item(self, item):
        self.items.remove(item)
        # self.calculate_total_price()

    def calculate_total_price(self):
        self.total_price = sum(item.price for item in self.items)

    @property
    def total_price(self):
        return sum(item.price for item in self.items)

    @total_price.setter
    def total_price(self, value):
        print('???')
        self.total_price = value


# 商品类
class Item:
    def __init__(self, name, price):
        self.name = name
        self.price = price


# 创建购物车
cart = Cart()

# 创建商品
item1 = Item("Item 1", 10)
item2 = Item("Item 2", 20)

# 将商品添加到购物车
cart.add_item(item1)
print("Total Price:", cart.total_price)  # 输出: Total Price: 10

cart.add_item(item2)
print("Total Price:", cart.total_price)  # 输出: Total Price: 30

# 从购物车中移除商品
cart.remove_item(item1)
print("Total Price:", cart.total_price)  # 输出: Total Price: 20

# 尝试使用setter，失败
# cart.total_price = 1
# print('total price2', cart.total_price)
