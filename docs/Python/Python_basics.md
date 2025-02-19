# Python 基础

## 基础语法

1. 字面量

字面量是指代码中直接写出的固定值，通常用于表示程序中的常量。在 Python 中，字面量可以通过 `print` 语句输出。

| **类型**       | **描述**                          | **说明**                                                     |
| -------------- | --------------------------------- | ------------------------------------------------------------ |
| **Number**     | `int`, `float`, `complex`, `bool` | 包括整数、浮点数、复数和布尔值（`True` 表示真，`False` 表示假）。 |
| **String**     | 描述文本的一种数据                | 字符串由任意数量的字符组成，用于表示文本。                   |
| **List**       | 有序的可变序列                    | Python 中使用最频繁的数据类型，用于有序记录一组数据。        |
| **Tuple**      | 有序的不可变序列                  | 可有序记录一组不可变的数据集合。                             |
| **Set**        | 无序不重复集合                    | 可无序记录一组不重复的数据。                                 |
| **Dictionary** | 无序键值对集合                    | 用于存储键值对（Key-Value）的无序数据集合。                  |

## 运算符

=== "算术运算符"

    Python 支持丰富的算术运算符，用于完成常见的数学计算。

    | **运算符** | **描述** | **实例**                                          |
    | ---------- | -------- | ------------------------------------------------- |
    | `+`        | 加       | 两个对象相加：`a + b`，例如 `10 + 20 = 30`        |
    | `-`        | 减       | 一个数减去另一个数：`a - b`，例如 `10 - 20 = -10` |
    | `*`        | 乘       | 两个数相乘：`a * b`，例如 `10 * 20 = 200`         |
    | `/`        | 除       | 返回商：`b / a`，例如 `20 / 10 = 2.0`             |
    | `//`       | 取整除   | 返回商的整数部分：`9 // 2 = 4`                    |
    | `%`        | 取余数   | 返回余数：`b % a`，例如 `20 % 10 = 0`             |
    | `**`       | 指数     | 计算幂：`a ** b`，例如 `10 ** 2 = 100`            |

=== "复合赋值运算符"

    复合赋值运算符是简化表达式的一种方式，将运算与赋值合并到一起。

    | **运算符** | **描述**         | **实例**                      |
    | ---------- | ---------------- | ----------------------------- |
    | `+=`       | 加法赋值运算符   | `c += a` 等效于 `c = c + a`   |
    | `-=`       | 减法赋值运算符   | `c -= a` 等效于 `c = c - a`   |
    | `*=`       | 乘法赋值运算符   | `c *= a` 等效于 `c = c * a`   |
    | `/=`       | 除法赋值运算符   | `c /= a` 等效于 `c = c / a`   |
    | `%=`       | 取余赋值运算符   | `c %= a` 等效于 `c = c % a`   |
    | `**=`      | 幂赋值运算符     | `c **= a` 等效于 `c = c ** a` |
    | `//=`      | 取整除赋值运算符 | `c //= a` 等效于 `c = c // a` |

## **4. 示例代码**

### 字面量

```python
# 数字字面量
a = 10          # 整数
b = 3.14        # 浮点数
c = True        # 布尔值

# 字符串字面量
text = "Hello, Python!"

# 列表字面量
my_list = [1, 2, 3]

# 字典字面量
my_dict = {"key1": "value1", "key2": "value2"}

# 输出各种字面量
print(a, b, c)
print(text)
print(my_list)
print(my_dict)
```

### 算术运算符

```python
a, b = 10, 20

# 加法
print(a + b)  # Output: 30

# 指数运算
print(a ** 2)  # Output: 100

# 取余和取整
print(b % a)   # Output: 0
print(b // a)  # Output: 2
```

### 复合赋值运算符

```python
c = 10

# 加法赋值
c += 5
print(c)  # Output: 15

# 幂赋值
c **= 2
print(c)  # Output: 225
```

## 5. 总结

- **字面量** 是程序中直接写出的固定值，用于构造变量或常量。
- **运算符** 提供了灵活的数学计算能力，复合赋值运算符可以简化代码。
- 通过熟练使用字面量和运算符，可以有效提升编程效率和代码可读性。

## 数据容器

数据容器是可以存储多个元素的 Python 数据类型。数据容器包括`list`（列表）`tuple`（元组）`str`（字符串）`set`（集合）`dict`（字典）。这些数据容器各有特点，但都具备存储多个元素的能力。

### 1. 列表（List）

```py
# 1.字面量
[元素1,元素2,元素3,元素4,元素5,...]
# 2.定义变量
变量名称 = [元素1,元素2,元素3,元素4,元素5,...]
# 3.定义空列表
变量名称 = []
变量名称 = list()
```

- 列表中的每个数据单元被称为 “元素”。

- 元素的数据类型没有限制，可以是任何类型，甚至可以是嵌套列表。

```py
# 定义一个列表
my_list = ["wuyang", "guojiaqi", "jiangjun"]
print(my_list)
print(type(my_list))

# 定义一个嵌套列表
my_list = [[1, 2, 3], [4, 5, 6]]
print(my_list)
```

### 索引

- **定义：** 列表中每个元素都有一个唯一的编号，称为索引。

- **索引规则：**

  - **正向索引：** 从前向后，从 `0` 开始递增。
  - **反向索引：** 从后向前，从 `-1` 开始递减。

  ```py
  # 通过索引取值
  列表[索引]
  ```

- **注意：** 索引超出范围会导致错误。

**示例：**

```py
my_list = ["wuyang", "guojiaqi", "jiangjun"]

# 正向索引
print(my_list[0])  # Output: wuyang
print(my_list[1])  # Output: guojiaqi

# 反向索引
print(my_list[-1])  # Output: jiangjun

# 嵌套列表索引
my_list2 = [[1, 2, 3], [4, 5, 6]]
print(my_list2[0])       # Output: [1, 2, 3]
print(my_list2[0][0])    # Output: 1
```

### 列表的方法

添加元素

- **`append`**：在列表尾部追加元素。
- **`extend`**：将另一个数据容器的元素追加到列表中。
- **`insert`**：在指定位置插入元素。

```py
my_list = [1, 2, 3]
my_list.append(4)  # [1, 2, 3, 4]
my_list.extend([5, 6])  # [1, 2, 3, 4, 5, 6]
my_list.insert(2, "new")  # [1, 2, 'new', 3, 4, 5, 6]
print(my_list)
```

#### **删除元素**

- **`del`**：通过索引删除元素。
- **`pop`**：删除并返回指定索引的元素（默认最后一个）。
- **`remove`**：删除列表中第一个匹配的值。
- **`clear`**：清空列表。

```py
my_list = [1, 2, 3, 2]
del my_list[1]  # [1, 3, 2]
removed = my_list.pop(0)  # 返回 1，列表变为 [3, 2]
my_list.remove(2)  # 删除第一个 2，结果：[3]
my_list.clear()  # []
```

#### **查找元素**

- **`index`**：返回元素的索引（若不存在抛异常）。
- **`count`**：统计元素在列表中出现的次数。

```py
my_list = ["wuyang", "guojiaqi", "jiangjun"]
print(my_list.index("guojiaqi"))  # 输出：1
print(my_list.count("jiangjun"))  # 输出：1
```

#### **统计与排序**

- **`len`**：返回列表长度。
- **`sort`**：对列表排序（默认升序，支持自定义规则）。
- **`reverse`**：反转列表顺序。

```py
my_list = [5, 2, 9, 1]
print(len(my_list))  # 输出：4
my_list.sort()  # [1, 2, 5, 9]
my_list.sort(reverse=True)  # [9, 5, 2, 1]
my_list.reverse()  # [1, 2, 5, 9]
```

###  列表的小技巧

#### 列表推导式

一种简洁创建列表的方法：

```py
# 创建一个平方列表
squares = [x**2 for x in range(5)]
print(squares)  # 输出：[0, 1, 4, 9, 16]

# 条件筛选
evens = [x for x in range(10) if x % 2 == 0]
print(evens)  # 输出：[0, 2, 4, 6, 8]
```

#### 解包操作

- 用 `*` 解包列表。

```py
numbers = [1, 2, 3, 4]
print(*numbers)  # 输出：1 2 3 4
```

## 函数

函数：是组织好的，可重复使用的，用来实现特定功能的代码段。

### 1. 函数的定义

```python
# 语法
def 函数名(传入参数)：
    函数体
    return 返回值

# 参数跟返回值不需要，可以省略。

# 示例
def greet(name):
    return f"Hello, {name}!"
# 调用
print(greet("ywu"))		# Output: Hello, ywu!
```

### 2. 函数的参数

位置参数

```py
# 定义函数
def add(x,y):
    result = x + y
    print(f"{x} + {y}的结果是：{result}"")

# 调用函数
add(5,6)	# 两个形参 x，y。 通过位置传递参数与形参一一对应。
```

默认参数

可变参数：`*args`

关键字参数：`**kwargs`

## Python的数据

### 变量相关

变量是Python中可以用来储存数值或者字符串的对象。变量有很多中类型，比如仅数值变量就细分有整数（int）、浮点数（float）等类型。而所有基于文本的变量都是字符串（str）类型。

#### 字符串

##### 字符串格式化

Python中的字符串格式化提供了多种方法来将变量插入到字符串中，以便生成易于阅读和使用的输出。以下是几种
常见的字符串格式化方法：

**1. 使用 `%` 操作符**

这是最早期的字符串格式化方法，类似于C语言中的`printf`。

```python
name = "John"
age = 30
formatted_string = "My name is %s and I am %d years old." % (name, age)
print(formatted_string)
```

有几个%?占位符，后面就跟几个变量或者值，顺序要对应好。如果只有一个%?，括号可以省略。

- 数字精度控制

 可以使用辅助符号“m .n”来控制数据的宽度和精度。

 	- %5d： 表示将整数宽度控制在5位。
 	- %5.2f：表示宽度控制为5，小数点精确2位。

##### 格式说明符

- `%s` 用于字符串
- `%d` 用于整数
- `%f` 用于浮点数
- `%%` 用于百分号

##### 2. 使用 `str.format()` 方法

这种方法引入了更强大的字符串格式化功能。

```python
name = "John"
age = 30
formatted_string = "My name is {} and I am {} years old.".format(name, age)
print(formatted_string)
```

你还可以通过位置和名称引用参数：

```python
formatted_string = "My name is {0} and I am {1} years old.".format(name, age)
print(formatted_string)

formatted_string = "My name is {name} and I am {age} years old.".format(name=name, age=age)
print(formatted_string)
```

##### 3. 使用 f-strings (格式化字符串字面量)

Python 3.6 引入了一种更简洁和高效的格式化方法：f-strings。

```python
name = "John"
age = 30
formatted_string = f"My name is {name} and I am {age} years old."
print(formatted_string)
```

你可以在 f-strings 中执行表达式：

```python
width = 10
height = 5
formatted_string = f"The area of the rectangle is {width * height} square units."
print(formatted_string)
```

##### 4. 使用 `string.Template` 类

这种方法提供了一种简单且安全的替代方案，适用于需要处理来自不可信源的数据的情况。

```python
from string import Template

template = Template("My name is $name and I am $age years old.")
formatted_string = template.safe_substitute(name="John", age=30)
print(formatted_string)
```

##### 详细格式化选项

在使用 `str.format()` 和 f-strings 时，可以使用格式说明符来控制字符串的显示。

##### 浮点数格式化

```python
number = 3.141592653589793
formatted_string = "{:.2f}".format(number)
print(formatted_string)

formatted_string = f"{number:.2f}"
print(formatted_string)
```

##### 对齐和填充

```python
text = "Hello"
formatted_string = "{:>10}".format(text)  # 右对齐
print(formatted_string)

formatted_string = "{:<10}".format(text)  # 左对齐
print(formatted_string)

formatted_string = "{:^10}".format(text)  # 居中对齐
print(formatted_string)

formatted_string = "{:*^10}".format(text)  # 居中对齐，并使用 * 填充
print(formatted_string)
```

##### 示例

以下是各种字符串格式化方法的综合示例：

```python
name = "Alice"
age = 28
balance = 1234.567

# 使用 % 操作符
print("Name: %s, Age: %d, Balance: %.2f" % (name, age, balance))

# 使用 str.format() 方法
print("Name: {}, Age: {}, Balance: {:.2f}".format(name, age, balance))

# 使用 f-strings
print(f"Name: {name}, Age: {age}, Balance: {balance:.2f}")

# 使用 string.Template
from string import Template
template = Template("Name: $name, Age: $age, Balance: $balance")
print(template.safe_substitute(name=name, age=age, balance=f"{balance:.2f}"))
```

这些方法各有优缺点，根据具体需求选择合适的方法可以提高代码的可读性和可维护性。



## 面向对象

判断对象的类型

```py
class A:
    pass


class B(A):
    pass


b = B()

# isinstance() 函数来判断一个对象是否是一个已知的类型，也参考类型继承关系
print(isinstance(b, B))
print(isinstance(b, A))

# == 是判断两个对象的值是否相等
print(type(b) == B)  # True
print(type(b) == A)  # False

# is 是判断两个对象的标识(内存地址)是否相等
print(type(b) is B)  # True
```



#### 魔术方法

`if __name__ == '__main__':` 是 Python 中一个常见的条件语句，它用于判断当前模块是否是直接被运行的。



