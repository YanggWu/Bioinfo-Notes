# 常见标准库

## OS

Python的`os`模块提供了一些与操作系统交互的功能，这使得程序能够执行如文件和目录的操作、获取系统信息、执行命令等。以下是一些常用的`os`模块函数的详细介绍：

### 文件和目录操作

1. `os.getcwd()`:
   - 获取当前工作目录的路径。
   - 示例：

```python
import os
current_directory = os.getcwd()
print(current_directory)
```

2. `os.chdir(path)`:
   - 改变当前工作目录。
   - 示例：

```python
os.chdir('/path/to/directory')
```

3. `os.listdir(path='.')`:
   - 返回指定目录下的文件和目录列表。
   - 示例：

```python
files = os.listdir('/path/to/directory')
print(files)
```

4. `os.mkdir(path)`:
   - 创建一个新的目录。
   - 示例：

```python
os.mkdir('/path/to/new_directory')
```

5. `os.makedirs(path)`:
   - 递归地创建目录。
   - 示例：

```python
os.makedirs('/path/to/new_directory/sub_directory')
```

6. `os.remove(path)`:
   - 删除文件。
   - 示例：

```python
os.remove('/path/to/file.txt')
```

7. `os.rmdir(path)`:
   - 删除目录。
   - 示例：

```python
os.rmdir('/path/to/directory')
```

8. `os.removedirs(path)`:
   - 递归地删除目录。
   - 示例：

```python
os.removedirs('/path/to/directory/sub_directory')
```

9. `os.rename(src, dst)`:
   - 重命名文件或目录。
   - 示例：

```python
os.rename('/path/to/file.txt', '/path/to/new_file.txt')
```

10. `os.stat(path)`:
   - 获取文件或目录的状态信息。
   - 示例：

```python
stats = os.stat('/path/to/file.txt')
print(stats)
```

1. `os.walk()`

   `os.walk` 是一个生成器函数，它返回一个三元组 `(dirpath, dirnames, filenames)`。具体来说：

   - `dirpath`：当前遍历到的目录的路径（字符串）。
   - `dirnames`：该目录下的子目录名列表（列表中的每一项是一个字符串）。
   - `filenames`：该目录下的文件名列表（列表中的每一项是一个字符串）。

### 路径操作

1. `os.path.basename(path)`:
   - 返回路径的基本名称。
   - 示例：

```python
base_name = os.path.basename('/path/to/file.txt')
print(base_name)  # 输出 'file.txt'
```

2. `os.path.dirname(path)`:
   - 返回路径的目录部分。
   - 示例：

```python
dir_name = os.path.dirname('/path/to/file.txt')
print(dir_name)  # 输出 '/path/to'
```

3. `os.path.join(path, *paths)`:
   - 将多个路径组合成一个路径。
   - 示例：

```python
full_path = os.path.join('/path', 'to', 'file.txt')
print(full_path)  # 输出 '/path/to/file.txt'
```

4. `os.path.split(path)`:
   - 将路径拆分为目录和文件名。
   - 示例：

```python
dir_name, base_name = os.path.split('/path/to/file.txt')
print(dir_name)   # 输出 '/path/to'
print(base_name)  # 输出 'file.txt'
```

5. `os.path.exists(path)`:
   - 检查路径是否存在。
   - 示例：

```python
exists = os.path.exists('/path/to/file.txt')
print(exists)  # 输出 True 或 False
```

6. `os.path.isabs(path)`:
   - 检查路径是否为绝对路径。
   - 示例：

```python
is_absolute = os.path.isabs('/path/to/file.txt')
print(is_absolute)  # 输出 True 或 False
```

7. `os.path.isfile(path)`:
   - 检查路径是否为文件。
   - 示例：

```python
is_file = os.path.isfile('/path/to/file.txt')
print(is_file)  # 输出 True 或 False
```

8. `os.path.isdir(path)`:
   - 检查路径是否为目录。
   - 示例：

```python
is_dir = os.path.isdir('/path/to/directory')
print(is_dir)  # 输出 True 或 False
```

9. `os.path.getsize(path)`:
   - 返回文件的大小（以字节为单位）。
   - 示例：

```python
size = os.path.getsize('/path/to/file.txt')
print(size)
```

10. `os.path.abspath(path)`:
   - 返回路径的绝对路径。
   - 示例：

```python
abs_path = os.path.abspath('file.txt')
print(abs_path)
```

### 系统信息

1. `os.name`:
   - 返回当前操作系统的名称。
   - 示例：

```python
os_name = os.name
print(os_name)  # 输出 'posix', 'nt' 等
```

2. `os.environ`:
   - 一个表示系统环境变量的字典。
   - 示例：

```python
home_dir = os.environ['HOME']
print(home_dir)
```

### 执行系统命令

1. `os.system(command)`:
   - 在子终端中执行系统命令。
   - 示例：

```python
os.system('ls -l')
```

2. `os.popen(command)`:
   - 执行系统命令并返回一个文件对象。
   - 示例：

```python
stream = os.popen('ls -l')
output = stream.read()
print(output)
```

`os`模块是Python标准库中功能非常强大的一个模块，它提供了对操作系统功能的访问，适用于需要进行文件操作、路径处理和执行系统命令的场景。

## pathlib

`pathlib` 是 Python 的标准库，提供了面向对象的方式进行路径操作，能够替代 `os.path` 进行文件和目录的管理，代码更简洁、可读性更强。

**创建路径对象：**

`Path` 是 `pathlib` 的核心类，用于表示和操作文件路径。

```py
from pathlib import Path

# 当前目录
p = Path(".")
print(p.resolve())  # 获取当前目录的绝对路径

# 指定路径
p = Path("/home/user/documents")
print(p)
```

### 获取路径信息

- `resolve()`: 返回路径的绝对路径，并去掉 `.` 和 `..`：

```py
p = Path(".")
print(p.resolve())  # /home/user/current_directory
```

- `name`：获取文件名（包括扩展名）

```py
p = Path("/home/user/file.txt")
print(p.name)  # file.txt
```

- `stem`：获取文件名（不包括扩展名）：

```py
print(p.stem)  # file
```

- `suffix`：获取文件扩展名：

```py
print(p.suffix)  # .txt
```

- `parent`：获取父目录：

```py
print(p.parent)  # /home/user
```



## sys

`sys`模块是Python标准库的一部分，提供了一些与Python解释器进行交互的函数和变量。这个模块常用于处理运行时环境配置、管理输入输出、操作Python运行时环境以及访问传递给脚本的命令行参数。

以下是`sys`模块的详细介绍：

### 1. **命令行参数**

- `sys.argv`: 这是一个列表，包含了运行Python脚本时传递的命令行参数。`sys.argv[0]`是脚本的名称，后面的元素是传递给脚本的参数。

```python
import sys

# 示例: python script.py arg1 arg2
print(sys.argv)  # 输出: ['script.py', 'arg1', 'arg2']
```

### 2. **程序退出与错误处理**

- `sys.exit([arg])`: 该函数用于退出程序，并可以选择性地返回一个状态码。通常，`0`表示程序正常退出，非零表示异常退出。你也可以传递一个字符串作为错误消息。

```python
import sys

# 正常退出
sys.exit(0)

# 异常退出
sys.exit("Error occurred")
```

### 3. **标准输入输出**

- `sys.stdin`: 表示标准输入（通常是键盘输入），你可以用它来从命令行或其他输入流中获取输入。

```python
import sys

# 从标准输入读取一行
user_input = sys.stdin.readline()
print(f"输入内容: {user_input}")
```

- `sys.stdout`: 表示标准输出（通常是控制台），可以用它来输出数据。

```python
import sys

# 打印到标准输出
sys.stdout.write("Hello, World!\n")
```

- `sys.stderr`: 表示标准错误输出流，通常用于输出错误信息。

```python
import sys

# 输出错误信息
sys.stderr.write("Error: Something went wrong\n")
```

### 4. **解释器相关**

- `sys.version`: 获取Python解释器的版本信息。

```python
import sys

print(sys.version)  # 输出Python版本信息
```

- `sys.platform`: 获取当前运行Python的操作系统平台名称。

```python
import sys

print(sys.platform)  # 输出: 'linux', 'darwin'（macOS）, 'win32' 等
```

- `sys.path`: 返回模块搜索路径的列表，Python会在这些路径中查找要导入的模块。可以动态地修改这个列表来改变模块的查找路径。

```python
import sys

print(sys.path)  # 输出当前的模块搜索路径列表
```

### 5. **内存管理**

- `sys.getsizeof(object)`: 返回对象在内存中的大小，单位是字节。

```python
import sys

x = [1, 2, 3, 4, 5]
print(sys.getsizeof(x))  # 输出对象x的字节大小
```

### 6. **运行时信息**

- `sys.modules`: 一个字典，包含了当前导入的所有模块的名称和对应的模块对象。

```python
import sys

print(sys.modules.keys())  # 输出当前已导入的模块名称
```

- `sys.maxsize`: 返回Python中整数的最大值（与系统的位数相关）。

```python
import sys

print(sys.maxsize)  # 输出系统支持的最大整数值
```

### 7. **异常处理**

- `sys.exc_info()`: 返回当前处理的异常类、异常实例以及异常回溯对象。这对于在异常处理过程中获取更多上下文信息非常有用。

```python
import sys

try:
    1 / 0
except:
    print(sys.exc_info())  # 输出异常信息
```

### 总结

`sys`模块为Python程序提供了与解释器、操作系统交互的多种方法，使得程序员能够更加灵活地控制程序的执行环境和行为。在编写复杂的Python脚本时，熟练使用`sys`模块会非常有帮助。