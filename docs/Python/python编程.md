## **装饰器**

在 Python 中，装饰器（decorator） 是一种用于修改或增强函数或类行为的高级特性。它本质上是一个函数，接受另一个函数或类作为输入，并返回一个经过修改或增强的新函数或类。

装饰器的核心功能是**在不修改原始代码的情况下，为函数或类添加新的功能**。以下是装饰器的详细介绍。

### 基本概念

**装饰器作用**：动态修改函数或类的行为，常用于日志记录、性能统计、权限控制、缓存等场景。

**核心结构**：

- 一个函数接受另一个函数作为输入。
- 返回一个新的函数。

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

**使用装饰器**

使用 `@` 语法将装饰器应用到函数或类。

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

