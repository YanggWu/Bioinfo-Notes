## **装饰器**

在 Python 中，装饰器（decorator） 是一种用于修改或增强函数或类行为的高级特性。它本质上是一个函数，接受另一个函数或类作为输入，并返回一个经过修改或增强的新函数或类。

装饰器的核心功能是在不修改原始代码的情况下，为函数或类添加新的功能。以下是装饰器的详细介绍。

### 基本概念

**装饰器作用**：动态修改函数或类的行为，常用于日志记录、性能统计、权限控制、缓存等场景。

**使用技巧**：

- 理解装饰器的嵌套与参数传递。
- 多实践典型场景（如日志、缓存、权限验证）。

```py
def decorator(func):
    def wrapper(*args, **kwargs):
        # 在调用原始函数之前的操作
        result = func(*args, **kwargs)  # 调用原始函数
        # 在调用原始函数之后的操作
        return result
    return wrapper
```

使用 `@` 语法将装饰器应用到函数或类：

```py
def simple_decorator(func):
    def wrapper():
        print("Before the function is called.")
        func()
        print("After the function is called.")
    return wrapper

@simple_decorator
def say_hello():
    print("Hello, World!")

say_hello()
```

### 装饰器的常见用途

1. **函数执行时间统计**

   ```py
   import time
   
   def timing_decorator(func):
       def wrapper(*args, **kwargs):
           start_time = time.time()
           result = func(*args, **kwargs)
           end_time = time.time()
           print(f"{func.__name__} executed in {end_time - start_time:.4f} seconds")
           return result
       return wrapper
   
   @timing_decorator
   def slow_function():
       time.sleep(2)
       print("Finished slow function")
   
   slow_function()
   ```

2. **权限验证**

   ```
   def require_permission(user):
       def decorator(func):
           def wrapper(*args, **kwargs):
               if user.get("is_admin"):
                   return func(*args, **kwargs)
               else:
                   print("Permission denied.")
           return wrapper
       return decorator
   
   current_user = {"name": "Alice", "is_admin": False}
   
   @require_permission(current_user)
   def delete_database():
       print("Database deleted!")
   
   delete_database()
   ```

3. **缓存结果**

   ```py
   def cache_decorator(func):
       cache = {}
       def wrapper(*args):
           if args in cache:
               print("Returning cached result")
               return cache[args]
           result = func(*args)
           cache[args] = result
           return result
       return wrapper
   
   @cache_decorator
   def add(a, b):
       return a + b
   
   print(add(3, 4))
   print(add(3, 4))  # 使用缓存
   ```

## 深入类和对象

## 抽象基类

抽象基类（Abstract Base Class，简称 ABC）是 Python 中的一种特殊类，用于定义抽象类（不能实例化的类），它可以作为其他类的基类。抽象基类用于强制要求子类实现某些特定方法，从而提供一种规范化的接口设计。抽象基类通过 `abc` 模块实现

**类中的方法：**

1. 实例方法（Instance Method）最常见的 Python 方法类型。第一个参数是self，它是定义在类中的普通方法，通常访问实例的属性和方法。
2. 类方法（Class Method）类方法是与类相关的方法，而不是与实例相关。它不需要实例对象就能调用。可以访问类的属性和方法，但不能直接访问实例的属性。第一个参数是 `cls`，表示类本身，通过 `@classmethod` 装饰器定义。
3. 静态方法（Static Method）静态方法不依赖于类的实例或类本身。它与类的实例无关，通常用于定义那些不需要访问类或实例属性的方法。静态方法没有 `self` 或 `cls` 参数，可以通过类或实例调用。

**私有属性：**

```py
# 私有属性
class Person:
    def __init__(self, name, money):
        self.name = name
        self.__money = money

    # 在类的内部可以访问私有属性
    def get_money(self):
        print(self.__money)


# 在外部无法使用私有属性
p = Person('张三', 100)
p.get_money()
```

**try 语句的使用**

```py
# try 语句的使用

try:
    print("程序运行。。。")
except KeyError:
    print("程序出现异常。。。")
else:
    print("程序未产生异常时则运行当前代码块。。。")
finally:
    print("程序无论是否出现异常都执行。。。")
```

函数中使用 try 语句注意 return返回的位置：如果finally中出现return，则优先返回finally中的值。

```py
# 在函数中使用try语句
def func_try_except():
    try:
        print("程序运行。。。")
        # raise KeyError
        return 1
    except KeyError:
        print("程序出现异常。。。")
        return 2
    else:
        print("程序未产生异常时则运行当前代码块。。。")
        return 3
    finally:
        print("程序无论是否出现异常都执行。。。")
        return 4


res = func_try_except()
print(res)
```

**上下文管理协议**

```py
# 上下文管理协议 - 魔术方法
class Sample:
    def __enter__(self):
        print("进入上下文管理器。。。")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        print("退出上下文管理器。。。")

    @staticmethod
    def run():
        print("程序运行。。。")


with Sample() as sample:
    sample.run()
```

