# Python基础

## 基础语法

### 字面量

> 代码中，被写在代码中的固定的值，称之为字面量。基于print语句完成各类字面量的输出。
>

**常用值的类型**

|    类型    | 描述                             | 说明                                                     |
| :--------: | :------------------------------- | -------------------------------------------------------- |
|   Number   | int    float      complex   bool | 整数，  浮点数，  复数， 布尔，True表示真，False表示假。 |
|   String   | 描述文本的一种数据               | 字符串（string）由任意数量的字符组成                     |
|    List    | 有序的可变序列                   | Python中使用最频繁的数据类型，可有序记录一堆数据         |
|   Tuple    | 有序的不可变序列                 | 可有序记录一堆不可变的Python数据集合                     |
|    Set     | 无序不重复集合                   | 可无序记录一堆不重复的Python数据集合                     |
| Dictionary | 无序Key-Value集合                | 可无序记录一堆Key-Value型的Python数据集合                |

### 算数运算符

| 运算符 | 描述   | 实例                                                         |
| :------: | :------: | ------------------------------------------------------------ |
| +      | 加     | 两个对象相加  a + b 输出结果  30                             |
| -      | 减     | 得到负数或是一个数减去另一个数  a - b 输出结果  -10          |
| *      | 乘     | 两个数相乘或是返回一个被重复若干次的字符串  a * b 输出结果  200 |
| /      | 除     | b /  a 输出结果 2                                            |
| //     | 取整除 | 返回商的整数部分  9//2 输出结果  4 , 9.0//2.0 输出结果  4.0  |
| %      | 取余   | 返回除法的余数  b % a 输出结果  0                            |
| **     | 指数   | a**b  为10的20次方，  输出结果 100000000000000000000         |

### 复合赋值运算符

| **运算符** | **描述**         | **实例**                    |
| :----------: | :------------- | --------------------------- |
| +=         | 加法赋值运算符   | c  += a 等效于 c =  c + a   |
| -=         | 减法赋值运算符   | c  -= a 等效于 c =  c - a   |
| *=         | 乘法赋值运算符   | c  *= a 等效于 c =  c * a   |
| /=         | 除法赋值运算符   | c  /= a 等效于 c =  c / a   |
| %=         | 取余赋值运算符   | c  %= a 等效于 c =  c % a   |
| **=        | 幂赋值运算符     | c  **= a 等效于 c =  c ** a |
| //=        | 取整除赋值运算符 | c  //= a 等效于 c =  c // a |

### 1. 字符串

- **字符串格式化**

  - 常用的三类占位符：
  
  | 格式符号 | 转化 |
  | :---: | :----------------------------- |
  |  %s  | 将内容转化为字符串，放入站位位置 |
  |  %d  |  将内容转化为整数，放入站位位置  |
  |  %f  | 将内容转化为浮点数，放入站位位置 |
  
	有几个%?占位符，后面就跟几个变量或者值，顺序要对应好。如果只有一个%?，括号可以省略。
	- 数字精度控制
	
	  可以使用辅助符号“m .n”来控制数据的宽度和精度。
	  
	  - %5d： 表示将整数宽度控制在5位。
	  
	  - %5.2f：表示宽度控制为5，小数点精确2位。

```py
print('Age: %s. Gender: %s' % (25, True))
print('growth rate: %d %%' % 7)
print('%2d-%02d' % (3, 1)) # %2d 表示多出两个空格，%02d 0个空格两位数数字
print('%.2f' % 3.1415926) # 保留两位小数
```

> 快速格式化可以使用
>
> f” {变量} {变量} “

- **表达式格式化**

  表达式：一条具有明确执行结果的代码语句

  <img src="https://raw.githubusercontent.com/wuayng-owl/Image/main/markdown_image/202312281828637.png" style="zoom:50%;" />

  在无需使用变量进行数据存储的时候，可以**直接格式化表达式**，简化代码。

```py
# 练习
name = "wuyang"
stock_price = 28.88
stock_code = "003456"
stock_price_daily_growth_factor = 1.2
growth_days = 7

print(f"公司：{name},股价代码：{stock_code},当前股价：{stock_price}")
print("每日增长系数是：%s,经过%d天的增涨后,股价达到了:%.2f" % 
      (stock_price_daily_growth_factor,
       growth_days,
       stock_price*growth_days*stock_price_daily_growth_factor))

#> 公司：wuyang,股价代码：003456,当前股价：28.88
#> 每日增长系数是：1.2,经过7天的增涨后,股价达到了:242.59
```

### 2. 数据输入 input

1. input()语句的功能是，获取键盘输入的数据
2. 可以使用：input(提示信息)，用以在使用者输入内容之前显示提示信息。
3. 要注意，无论键盘输入什么类型的数据，获取到的数据**==永远都是字符串类型==**

```py
#练习
v1=input("请输入一个字符串：")
v2=input("请输入一个整数：")
v3=input("请输入一个浮点数：")
v4=input("请输入一个布尔类型：")
print(f"{type(v1)},内容是：{v1}")
print(f"{type(v2)},内容是：{v2}")
print(f"{type(v3)},内容是：{v3}")
print(f"{type(v4)},内容是：{v4}")

#> <class 'str'>,内容是：wuyang
#> <class 'str'>,内容是：5
#> <class 'str'>,内容是：5.4
#> <class 'str'>,内容是：Fasle
```

## Python判断语句

### 1. 布尔类型和比较运算符

布尔类型的字面量只有：True，False

可以定义变量存储布尔类型数据，也可通过比较运算符得到。

|  **运算符**  |                         **描述**                           |            **示例**             |
| :----------------: | :------------------------------------------------- | :---------------------------: |
|     ==     |         判断内容是否相等，满足为True，不满足为False          | 如a=3,b=3，则(a  == b)  为  True |
|     !=     |        判断内容是否不相等，满足为True，不满足为False         | 如a=1,b=3，则(a  != b) 为  True  |
|     >      |  判断运算符左侧内容是否大于右侧  满足为True，不满足为False   |  如a=7,b=3，则(a  > b)  为 True  |
|     <      |  判断运算符左侧内容是否小于右侧  满足为True，不满足为False   |  如a=3,b=7，则(a  < b)  为 True  |
|     >=     | 判断运算符左侧内容是否大于等于右侧  满足为True，不满足为False | 如a=3,b=3，则(a  >= b) 为  True  |
|     <=     | 判断运算符左侧内容是否小于等于右侧  满足为True，不满足为False | 如a=3,b=3，则(a  <= b) 为  True  |

### 2. if语句

1. if语句的基本格式

```py
if 要判断的条件 :
	条件成立时，要做的事情
```

2. if语句的注意事项：

   1. 判断条件的结果一定要是布尔类型

   2. 不要忘记判断条件后的： 冒号

   3. 归属于if语句的代码块，需在前方填充4个空格缩进

```py
# 练习
print("欢迎来到儿童游乐园，儿童免费，成人收费。")
age = input("请输入你的年龄：")
if int(age) >= 18 :
    print("您已成年，游玩需要补票10元\n祝您游玩愉快")

#> 欢迎来到儿童游乐园，儿童免费，成人收费。
#> 请输入你的年龄：25
#> 您已成年，游玩需要补票10元
#> 祝您游玩愉快
```

### 3. if else语句

```py
if 要判断的条件 :
	条件成立时，要做的事情1
	条件成立时，要做的事情2
	...省略...
else:
	不满足条件时要做的事情1
	不满足条件时要做的事情2
	...省略...
```

1. if else 语句，其中

   - if和其代码块，条件满足时执行

   - else搭配if的判断条件，当不满足的时候执行

2. if else语句的注意事项：

   - else不需要判断条件，当if的条件不满足时，else执行

   - else的代码块，同样要4个空格作为缩进

```py
# 练习
print("欢迎来到动物园")
height = int(input("请输入你的身高（cm）："))
if height >= 120:
    print("您的身高超过120cm，游玩需要购票10元。")
else:
    print("您的身高未超出120cm，可以免费游玩。")

print("祝您游玩愉快。")
```

### 4. if elif else语句

```py
if 条件1:
	条件1满足应该做的事情1
	条件1满足应该做的事情2
	...省略...
elif 条件2:
	条件2满足应该做的事情1
	条件2满足应该做的事情2
	...省略...
elif 条件N:
	条件N满足应该做的事情1
	条件N满足应该做的事情2
	...省略...
else:
	所有条件都不满足应该做的事情
	......
```

1. if elif else语句的作用是？
   - 可以完成多个条件的判断

2. 使用if elif else的注意点有：

   - elif可以写多个

   - 判断是互斥且有序的，上一个满足后面的就不会判断了

   - 可以在条件判断中，直接写input语句，节省代码量

```py
# 练习
ideal_num = 888
if int(input("请输入第一次猜想的数字：")) == ideal_num:
    print("恭喜你第一次就猜对了")
elif int(input("不对，在猜一次：")) == ideal_num:
    print("恭喜你猜对了")
elif int(input("不对，在猜最后一次：")) == ideal_num:
    print("恭喜你最后一次猜对了")
else:
    print(f"Sorry,全部猜错了，我想的是：{ideal_num}")
```

### 5. 判断语句的嵌套

基础语法格式

```py
if 条件1:
	满足条件1 做事情
	
	if 条件2:
		满足条件2，做的事情
```

1. 嵌套判断语句可以用于多条件、多层次的逻辑判断

2. 嵌套判断语句可以根据需求，自由组合if elif else来构建多层次判断

3. 嵌套判断语句，一定要注意空格缩进，Python通过空格缩进来决定层次关系

### 实战案例

要求：

1. 数字随机产生，范围1-10

2. 有3次机会猜测数字，通过3层嵌套判断实现
3. 每次猜不中，会提示大了或小了

```py
import random
num = random.randint(1, 10)

num_1 = int(input("请输入第一次猜想的数字："))

if num_1 == num:
    print("恭喜你第一次就猜对了")
else:
    if num_1 > num:
        print("你猜的数字大了")
    else:
        print("你猜小了。")
    num_2 = int(input("再次输入你要猜测的数字："))

    if num_2 == num:
        print("恭喜你第二次猜中了")
    else:
        if num_2 > num:
            print("你猜的数字大了")
        else:
            print("你猜小了。")

        num_3=int(input("第三次输入你要猜测的数字："))
        if num_3 == num:
            print("恭喜第三次猜对了")
        else:
            print("三次机会都没有猜中！")
```

## Python循环语句

### 1. while循环语句

**while循环的基础**

```py
while 条件：
	条件满足时，做的事情1
	条件满足时，做的事情2
	条件满足时，做的事情3
	.....
# 只要条件满足，会无限循环执行。
```

1. while的条件需得到布尔类型，True表示继续循环，False表示结束循环
2. 需要**设置循环终止的条件**，如i += 1配合 i < 100，就能确保100次后停止，否则将无限循环
3. 空格缩进和if判断一样，都需要设置

```py
# while cycle practice
i = 0
sum = 0
while i < 100:
    i += 1
    sum += i
print(f"1-100累加的最终结果是:{sum}")
```

练习案例

```py
# while cycle practice 2
import random
num = random.randint(1,100)
i = 1
g_num = int(input("请输入你猜测的数字："))
while g_num != num:
    if g_num > num:
        print(f"第{i}次猜错了，你猜大了")
    else:
        print(f"第{i}次猜错了，你猜小了")
    g_num = int(input("请重新输入你猜的数字："))
    i += 1
print(f"恭喜你第{i}次猜对了数字：{num}")
```

**while循环的嵌套**

```py
while 条件1:
	条件1满足时，做的事情
	......
	
	while 条件2:
		条件2满足时，做的事情
		.....
```

案例练习

打印九九乘法表

```py
# while cycle practice
i = 1
while i <= 9:
    j = 1
    while j <= i:
        print(f"{j}*{i}={j * i}\t", end='')
        j += 1
    print("")
    i += 1
```

### 2. for循环

基础语法

```py
for 临时变量 in 待处理数据集:
	执行的代码
```

待处理数据集，严格来说，称之为：序列类型，指其内容可以一个个依次取出的一种类型。

包括：字符串、列表、元组、字典等

**range语句**

```py
rang(num) # 获取一个从0开始，到num结束的数字序列（不含num本身）
rang(num1,num2) # 获取一个从0开始，到num结束的数字序列（不含num本身）
rang(num1,num2,step)
```

### 案例练习

while

```py
# while cycle practice
import random
num = random.randint(1,100)
i = 1
g_num = int(input("请输入你猜测的数字："))
while g_num != num:
    if g_num > num:
        print(f"第{i}次猜错了，你猜大了")
    else:
        print(f"第{i}次猜错了，你猜小了")
    g_num = int(input("请重新输入你猜的数字："))
    i += 1
print(f"恭喜你第{i}次猜对了数字：{num}")
```

for

```py
# for cycle practice
seqs = "ATTTAAGGCCGTAGCGATTGGATGATTCGAT"
i = 0
for seq in seqs:
    if seq == "A":
        i += 1
print(f"这个序列中共含有：{i}个A")

# Combined training
money = 10000
for i in range(1,21):
    import random
    num = random.randint(1,10)
    if money <= 0:
        print("工资发完了，等下个月吧")
        break
    if num < 5:
        print(f"员工{i},绩效分{num}，低于5，不发工资，下一位")
        continue
    elif money > 0:
        money -=1000
        print(f"员工{i}满足，发放工资1000元，账户余额还剩余{money}")
```



## 函数

> 函数：是组织好的，可重复使用的，用来实现特定功能的代码段。

### 1. 函数的定义

```py
# 语法
def 函数名(传入参数)：
	函数体
	return 返回值

# 参数跟返回值不需要，可以省略。
```

### 2. 函数的参数

```py
# 定义函数
def add(x,y):
	result = x + y
	print(f"{x} + {y}的结果是：{result}"")
	
# 调用函数
add(5,6)
```

练习案例

```py
# def practice
def check(x):
    print("欢迎来到黑马程序员！请出示健康码，并配合测量体温.")
    if x <= 37.5:
        print(f"您的体温是：{x}度，体温正常请进！")
    else:
        print(f"您的体温是：{x}度，需要隔离")

check(39)
check(36)
```

### 3. 函数的返回值

1. 语法格式

```py
def 函数(参数...):
	函数体
	return 返回值

变量 = 函数(参数)

# practice
def add(x,y):
    result = x + y
    return result

r = add(10,26)
print(r)
```

2. None类型

   Python中有一个特殊的字面量：None，其类型是：<class 'NoneType'>

   无返回值的函数，实际上就是返回了：None这个字面量

   None表示：空的、无实际意义的意思

   函数返回的None，就表示，这个函数没有返回什么有意义的内容。

   也就是返回了空的意思。

### 4. 函数说明文档

语法格式

```py
def add(x,y):
    """
    add函数可以接收两个参数，进行2数相加的功能
    :param x:形参x表示其中一个数字
    :param y:形参y表示另一个数字
    :return:返回两数相加的结果
    """
    result = x + y
    print(f"2数相加的结果是：{result}")
    return result

r = add(5,6)
```

### 5. 函数的嵌套

所谓函数嵌套调用指的是一个函数里面又调用了另外一个函数

```py
def func_b():
	print("---2---")
	
def func_a():
	print("---1---")
	
	func_b
	
	print("---3---")
# 调用函数
func_a()
```

### 6. 变量的作用域

> 变量作用域指的是变量的作用范围（变量在哪里可用，在哪里不可用）
>
> 主要分为两类：局部变量和全局变量

1. **局部变量**

   所谓局部变量是定义在函数体内部的变量，即只在函数体内部生效

   ```py
   def test_a():
       num = 100
       print(num)
   
   test_a() # 100
   print(num) # 报错：name 'num' is not defined
   ```

   局部变量的作用：在函数体内部，临时保存数据，即当函数调用完成后，则销毁局部变量

2. **全局变量**

   所谓全局变量，指的是在函数体内、外都能生效的变量

   ```py
   # 定义全局变量
   num = 100
   
   def test_a():
       print(num) # 访问全局变量num，并打印变量num储存的数据
       
   test_a()
   ```

3. **global关键字**

   使用 global关键字 可以在函数内部声明变量为全局变量

   ```py
   num = 100
   def test_a():
       print(num)
   
   def test_b():
       # global 关键字声明num 是全局变量
       global num
       num = 200
       print(num)
   
   test_a()
   test_b()
   print(f"全局变量num = {num}")
   ```

### 7. 综合案例

```py
# combined training
money = 5000000
name = input("请输入您的姓名:")
i = 0
# 定义查询余额函数
def inquire():
    print(f"------------查询余额------------\n{name},您好，您的账户余额为：{money}元")
# 存款函数
def save(x):
    global money
    money += x
    print(f"------------存款-----------\n{name},您好,存款{x}元成功\n当前余额{money}")
# 取款函数
def withdraw(x):
    global money
    money -= x
    print(f"------------取款-----------\n{name},您好取款{x}元成功\n当前余额{money}")
# 主菜单函数
def menu():
    print(f"------------主菜单-----------\n{name},您好，欢迎来到银行ATM，请选择操作：")
    print("查询余额\t【输入1】\n存款\t\t【输入2】\n取款\t\t【输入3】\n退出\t\t【输入4】")
    global i
    i = int(input("请输入您的选择"))
    if i == 1:
        inquire()
    elif i == 2:
        save(int(input("请输入存款金额：")))
    elif i == 3:
        withdraw(int(input("请输入取款金额：")))
    else:
        i = 4
        print("您好，当前已退出主菜单")

while i != 4:
    menu()
```

## 数据容器

> 1. 什么是数据容器？
>
>    一种可以存储多个元素的Python数据类型
>
> 2. Python有哪些数据容器？
>
>    list(列表)、tuple(元组)、str(字符串)、set(集合)、dict(字典)
>
>    它们各有特点，但都满足可容纳多个元素的特性。

### 1. 列表

#### 1.1**基本语法**

```py
# 字面量
[元素1,元素2,元素3,元素4,元素5,...]
# 定义变量
变量名称 = [元素1,元素2,元素3,元素4,元素5,...]
# 定义空列表
变量名称 = []
变量名称 = list()
```

​	数据容器内的每一份数据，都称之为元素

​	元素的数据类型没有任何限制，甚至元素也可以是列表，这样就定义了嵌套列表

```py
# 定义一个列表
my_list = ["wuyang","guojiaqi","jiangjun"]
print(my_list)
print(type(my_list))
# 定义一个嵌套列表
my_list = [ [1,2,3],[5,6,7] ]
```

#### 1.2 **列表的索引**

> 1. 列表的下表索引是什么
>
>    列表的每一个元素，都有编号称之为下标索引
>
>    从前向后的方向，编号从0开始递增
>    
>    从后向前的方向，编号从-1开始递减
>    
> 2. 如何通过下标索引取出对应位置的元素呢？
>
>    列表[下标]，即可取出
>
> 3. 下标索引的注意事项：
>
>    要注意下标索引的取值范围，超出范围无法取出元素，并且会报错
>

```py
my_list = ["wuyang", "guojiaqi", "jiangjun"]
# 从前向后0开始
print(my_list[0])
print(my_list[1])
# 倒序取出--从-1开始
print(my_list[-1])

# 取出嵌套列表元素
my_list2 = [[1, 2, 3], [4, 5, 6]]
print(my_list2[0])
print(my_list2[0][0])
```

#### 1.3 **列表的查询功能（方法）**

> 函数是一个封装的代码单元，可以提供特定功能。
>
> 在Python中，如果将函数定义为class（类）的成员，那么函数会称之为：方法

**index方法**

```py
# 语法
列表.index(元素)

my_list = ["wuyang", "guojiaqi", "jiangjun"]
print(my_list.index("wuyang"))  # 结果：1
```

**列表修改**

```py
# 修改特定位置（索引）的元素值：
列表[下标] = 值

# 正向下标
mylist = [1, 2, 3]
mylist[0] = 5
print(mylist)
# 反向下标
mylist = [1, 2, 3]
mylist[-3] = 5
print(mylist)
```

**插入元素**

```py
列表.insert(下标, 元素) #在指定的下标位置，插入指定的元素

mylist = [1, 2, 3]
mylist.insert(0,"wuyang")
print(mylist)
```

**追加元素**

```py
# 将指定元素，追加到列表的尾部
列表.append(元素)

# 追加元素方式2，将其他数据容器的内容取出，依次追加
列表.extend(其它数据容器)

mylist = ["wuyang", "guojiaqi", "jiangjun", 666]
# 追加新元素到列表
print(f"追加元素后的列表：{mylist}")
# 追加列表
mylist = ["wuyang", "guojiaqi", "jiangjun"]
mylist2 = [1, 2, 3]
mylist.append(mylist2)
print(f"追加列表后的结果：{mylist}")
# 追加新列表中的元素
mylist = ["wuyang", "guojiaqi", "jiangjun"]
mylist2 = [1, 2, 3]
mylist.extend(mylist2)
print(f"追加列表后的结果：{mylist}")
```

**删除元素**

1. 通过下标删除元素

   ```py
   # 通过下标删除元素
   mylist = [1, 2, 3]
   del mylist[0]
   print(mylist)
   # 通过下标删除，并取出。
   mylist = [1, 2, 3]
   num = mylist.pop(0)
   print(f"取出后：{mylist},取出的元素是：{num}")
   ```

2. 直接删除元素

   ```py
   # 删除某元素在列表中的第一个匹配项
   列表.remove(元素)
   
   mylist = ["wuyang", "guojiaqi", "guojiaqi", "jiangjun"]
   mylist.remove("guojiaqi")
   print(mylist)
   ```

3.  清空列表

   ```py
   # 清空列表内容
   列表.clear()
   
   mylist = ["wuyang", "guojiaqi", "guojiaqi", "jiangjun"]
   mylist.clear()
   print(mylist)
   ```

**统计元素**

1. 统计某元素的数量

   ```py
   列表.count(元素)
   
   mylist = ["wuyang", "guojiaqi", "guojiaqi", "jiangjun"]
   count = mylist.count("guojiaqi")
   print(count)
   ```

2. 统计列表有多少元素

   ```py
   len(列表)
   
   mylist = ["wuyang", "guojiaqi", "guojiaqi", "jiangjun"]
   count = len(mylist)
   print(count)
   ```

#### 1.4 **练习案例**

```py
# list pracitce

mylist = [21, 25, 21, 23, 22, 20]
# 追加元素
mylist.append(31)
print(mylist)
# 追加新列表中的元素
mylist2 = [29, 33, 30]
mylist.extend(mylist2)
print(mylist)
# 取出第一个元素
mylist[0]
print(mylist[0])
# 取出最后一个元素
mylist[-1]
print(mylist[-1])
# 查找元素31位置
mylist.index(31)
print(mylist.index(31))
```

#### 1.5 列表方法总览

| **编号** | **使用方式**            | **作用**                                       |
| -------- | ----------------------- | ---------------------------------------------- |
| 1        | 列表.append(元素)       | 向列表中追加一个元素                           |
| 2        | 列表.extend(容器)       | 将数据容器的内容依次取出，追加到列表尾部       |
| 3        | 列表.insert(下标, 元素) | 在指定下标处，插入指定的元素                   |
| 4        | del 列表[下标]          | 删除列表指定下标元素                           |
| 5        | 列表.pop(下标)          | 删除列表指定下标元素                           |
| 6        | 列表.remove(元素)       | 从前向后，删除此元素第一个匹配项               |
| 7        | 列表.clear()            | 清空列表                                       |
| 8        | 列表.count(元素)        | 统计此元素在列表中出现的次数                   |
| 9        | 列表.index(元素)        | 查找指定元素在列表的下标  找不到报错ValueError |
| 10       | len(列表)               | 统计容器内有多少元素                           |

#### 1.6 列表的遍历

> 1. 什么是遍历？
>
> 将容器内的元素依次取出，并处理，称之为遍历操作
>
> 2. 如何遍历列表的元素？
>
> 可以使用while或for循环

1. **while循环**

   ```py
   mylist = ["wuyang", "guojiaqi", "jiangjun"]
   index = 0
   while index < len(mylist):
       element = mylist[index]
       print(f"列表中的元素有：{element}")
       index += 1
   ```

2. **for循环**

   ```py
   mylist = ["wuyang", "guojiaqi", "jiangjun"]
   for index in mylist:
       print(f"列表中的元素有：{index}")
   ```

3. **练习案例**

   ```py
   # combined training
   
   # 1.for cycle
   mylist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
   newlist = []
   for i in mylist:
       num = i % 2
       if num == 0:
           newlist.append(i)
   print(f"查看新列表：{newlist}")
   
   # 2.while cycle
   mylist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
   newlist = []
   index = 0
   while index < len(mylist):
       element = mylist[index]
       num = element % 2
       if num == 0:
           newlist.append(element)
       index += 1
   print(f"查看新列表：{newlist}")
   ```


### 2. 元组

> 元组同列表一样，都是可以封装多个、不同类型的元素在内。
>
> 但最大的不同点在于：
>
> **元组一旦定义完成，就不可修改**

```py
# 定义元组字面量
(元素1，元素2，....)
# 定义元组变量
变量名称 = (元素1,元素2,元素2	...)
# 定义空元组
变量 = ()			#方式1
变量 = tuple()	#方式2

# 注意事项
# 定义3个元素的元组
t1 = (1,"Hello",True)
# 定义1个元素的元组
t2 = ("Hello",) # 注意，必须带有逗号，否则不是元组类型
```

1. **元组的相关操作**

   元组由于不可修改的特性，操作方法非常少。

   | **编号** | **方法**  | **作用**                                           |
   | -------- | --------- | -------------------------------------------------- |
   | 1        | index()   | 查找某个数据，如果数据存在返回对应的下标，否则报错 |
   | 2        | count()   | 统计某个数据在当前元组出现的次数                   |
   | 3        | len(元组) | 统计元组内的元素个数                               |

   ```py
   # 不可修改元组的内容，否则会直接报错。
   t1 = (1, 2, 3)
   t1[0] = 5  # 报错
   
   # 但是可以修改元组内的list的内容（修改元素，增加，删除，反转等）
   ```

2. 练习案例

   ```py
   # tuple practice
   t1 = ('周杰伦', 11, ['football', 'music'])
   # 查询下标位置
   index = t1.index(11)
   print(index)
   # 删除学生姓名
   print(f'该学生的姓名是：{t1[0]}')
   # 删除爱好中的football
   del t1[2][0]
   print(f"删除后的元组：{t1}")
   # 增加爱好 coding 到list内
   t1[2].append("coding")
   print(f"增加后的元组：{t1}")
   ```

### 3. 字符串

> 尽管字符串看起来并不像：列表、元组那样，一看就是存放了许多数据的容器。
>
> 但不可否认的是，字符串同样也是数据容器的一员。
>
> 字符串是字符的容器，一个字符串可以存放任意数量的字符。

1. 字符串的下标（索引）

   和其它容器如：列表、元组一样，字符串也可以通过下标进行访问

   •从前向后，下标从0开始

   •从后向前，下标从-1开始

   ```py
   # 通过下标获取特定位置字符
   my_str = "wuyang"
   print(my_str[0])
   print(my_str[-1])
   # 同元组一样，字符串是一个：无法修改的数据容器
   ```

2. 查找特定字符串的下标索引值

   语法：字符串.index(字符串)

   ```py
   my_str = "wuyang and guojiaqi"
   print(my_str.index("and"))
   ```

3. 字符串的替换

   语法：字符串.replace(字符串1，字符串2）

    功能：将字符串内的全部：字符串1，替换为字符串2

    注意：不是修改字符串本身，而是得到了一个新字符串

   ```py
   mystr = "learn Linux"
   newstr = mystr.replace("Linux","python")
   print(newstr)
   # 可以看到,字符串mystr本身并没有发生变化
   ```

4. 字符串的分割

   语法：字符串.split(分隔符字符串）

   功能：按照指定的分隔符字符串，将字符串划分为多个字符串，并存入列表对象中

   注意：**字符串本身不变，而是得到了一个列表对象**

   ```py
   # character split
   my_str = "learning python bioinformatics"
   my_lst = my_str.split(" ")
   print(f"将字符串{my_str}切割后,得到一个列表对象:{my_lst}")
   ```

5. 字符串的规整操作

   语法:字符串.strip(字符串)

   ```py
   # .strip() 默认去除前后空格
   my_str = "  learning python bioinformatics      "
   new_str = my_str.strip(" ")
   print(f"取出前后空格得到的新字符串是:{new_str}")
   
   # .strip(字符串)
   my_str = "12learning python bioinformatics21"
   new_str = my_str.strip("12 ")
   print(f"新字符串是:{new_str}")
   # 传入的是“12” 其实就是：”1”和”2”都会移除，是按照单个字符
   ```

6. 统计字符串长度、某字符出现次数

   语法: len(字符串)

   语法：字符串.count(字符串)

   ```py
   my_str = "learning python bioinformatics"
   print(len(my_str))
   
   # .count() practice
   my_seq = "AAATTTGGCCAGCTTTGCAAA"
   num = my_seq.count("A")
   print(num)
   ```

7. 字符串常用操作

| **编号** | **操作**                             | **说明**                                                     |
| -------- | ------------------------------------ | ------------------------------------------------------------ |
| 1        | 字符串[下标]                         | 根据下标索引取出特定位置字符                                 |
| 2        | 字符串.index(字符串）                | 查找给定字符的第一个匹配项的下标                             |
| 3        | 字符串.replace(字符串1, 字符串2)     | 将字符串内的全部字符串1，替换为字符串2   不会修改原字符串，而是得到一个新的 |
| 4        | 字符串.split(字符串)                 | 按照给定字符串，对字符串进行分隔  不会修改原字符串，而是得到一个新的列表 |
| 5        | 字符串.strip()  字符串.strip(字符串) | 移除首尾的空格和换行符或指定字符串                           |
| 6        | 字符串.count(字符串)                 | 统计字符串内某字符串的出现次数                               |
| 7        | len(字符串)                          | 统计字符串的字符个数                                         |

> 字符串特点:
>
> •只可以存储字符串
>
> •长度任意（取决于内存大小）
>
> •支持下标索引
>
> •允许重复字符串存在
>
> •不可以修改（增加或删除元素等）
>
> •支持for循环

8. 练习案例

   ```py
   # character practice
   my_str = "itheima itcast boxuegu"
   num = my_str.count("it")
   new_str = my_str.replace(" ", "|")
   lst_str = new_str.split("|")
   print(f"字符串{my_str}中有: {num}个it字符")
   print(f"字符串{my_str}被替换空格后,结果:{new_str}")
   print(f"字符串{new_str},按照|分割后,得到{lst_str}")
   ```

### 4. 数据容库的切片

1. 序列

   序列是指：内容连续、有序，可使用下标索引的一类数据容器

   列表、元组、字符串，均可以可以视为序列。

2. 序列的切片 

   ```py
   # 语法
   序列[起始下标:结束下标:步长]
   
   # practice
   my_str = "万过薪月，员序程马黑来，nohtyP学"
   new_str = my_str[9:4:-1]
   print(new_str)
   new_lst = my_str.split("，")
   print(new_lst[1].replace("来","")[::-1])
   ```

### 5. 集合

> 集合元素无序、不可重复
>
> 所以集合不支持下标索引

```py
# 基本语法
# 定义集合字面量
{元素,元素,元素,....}
# 定义集合变量
变量名称 = {元素,元素,元素,....}
# 定义空集合
变量名称 = set()
```

#### 集合常用操作

```
# 添加新元素
集合.add(element)

# 移除元素
集合.remove(element)

# 从集合中随机
```



### 6. 字典

```py
# 定义字典
my_dict = {"wuyang": 99, "guojiaqi": 88, "jiangjun": 77}
# 定义空字典
my_dict2 = {}
my_dict3 = dict()
print(f"字典1的内容是：{my_dict}:类型：{type(my_dict)}")
print(f"字典2的内容是：{my_dict2}:类型：{type(my_dict2)}")

# 从字典的中基于key获取Value
my_dict = {"wuyang": 99, "guojiaqi": 88, "jiangjun": 77}
value = my_dict["wuyang"]
print(f"wuyang的值是：{value}")
```

#### 字典的嵌套

```
# 定义嵌套字典

```

## 函数

### 函数作为参数传递

```py
def test_func(compute):
    result = compute(1, 2)  # 确定compute是函数
    print(f"compute参数的类型是:{type(compute)}")
    print(f"计算的结果是:{result}")

# 定义一个函数,准备作为参数传入另一个函数
def compute(x, y):
    return x + y
# 调用,并传入函数
test_func(compute)
```

### lambda匿名函数

## 文件操作

### 文件读取

```py
# 打开文件
f = open("/Users/wuyang/Desktop/test_norm.txt", "r",encoding='UTF-8')
print(type(f))
# 读取文件 read()
print(f"读取10个字节的结果：{f.read(10)}")
print(f"read方法读取全部内容的结果是{f.read()}")

# 读取文件 readLines()
lines = f.readlines()  # 读取文件的全部行，封装到列表中
print(f"lines对象的类型是：{type(lines)}")
print(f"lines对象的内容是：{lines}")

# readline() # 一次读取一行内容
line1 = f.readline()
line2 = f.readline()
line3 = f.readline()
print(f"第一行数据是：{line1}")
print(f"第二行数据是：{line2}")
print(f"第三行数据是：{line3}")

# for 循环读取文件行
for line in f:
    print(f"每一行数据是：{line}")

# 文件关闭 .close() 解除对文件的占用
f.close()
```

