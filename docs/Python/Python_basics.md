### Python基础

## Python的数据结      

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

#### 魔术方法

`if __name__ == '__main__':` 是 Python 中一个常见的条件语句，它用于判断当前模块是否是直接被运行的。让我们来详细解释一下它的作用和原理：

##### `__name__` 变量

在 Python 中，`__name__` 是一个特殊的变量，用来表示当前模块的名字。当 Python 解释器运行一个脚本时，它会为每个模块设置一个特定的 `__name__` 变量值：

- 如果一个模块是直接被运行的，`__name__` 的值会是 `'__main__'`。
- 如果一个模块是被其他模块导入的，`__name__` 的值会是该模块的名字。

##### `if __name__ == '__main__':` 的作用

`if __name__ == '__main__':` 的主要作用是用来判断当前模块是否作为主程序运行。这样的条件语句通常用于以下两种情况：

1. **模块作为主程序运行时的逻辑处理**：
 - 当一个 Python 文件被直接运行时，解释器会把 `__name__` 设置为 `'__main__'`，然后执行 `if __name__ == '__main__':` 之后的代码块。
 - 这样可以在模块作为主程序运行时执行一些特定的初始化、逻辑处理或者测试代码。
2. **模块被其他模块导入时避免不必要的代码执行**：
 - 如果一个模块被其他模块导入，Python 解释器会将 `__name__` 设置为该模块的名字，而不是 `'__main__'`。
 - 使用 `if __name__ == '__main__':` 条件语句可以避免在被导入时执行不需要的初始化代码，因为只有当模块作为主程序运行时，才会执行 `if __name__ == '__main__':` 下的代码块。

##### 示例说明

 python
 # 模块中的一些定义和函数

 def main():
  # 主程序逻辑
  print("This is the main function")

 # 判断模块是否作为主程序运行
 if __name__ == '__main__':
  main()

在这个示例中：

- `main()` 函数定义了主程序的逻辑处理。
- `if __name__ == '__main__':` 判断当前模块是否作为主程序运行。
- 如果该模块作为主程序运行，那么调用 `main()` 函数执行主程序的逻辑。

通过这种方式，可以使 Python 脚本既可以作为独立的程序执行，又可以作为模块被其他程序导入使用，保证了代码的灵活性和可重用性。
