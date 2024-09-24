# 正则表达式

当使用 Python 中的 `re` 标准库时，通常用于处理正则表达式。正则表达式是一种强大的文本匹配和搜索工具，它允许你通过模式描述来查找、替换和操作字符串。

[正则表达式指南](https://docs.python.org/zh-cn/3/howto/regex.html#regex-howto)

| **方法 / 属性** | **目的**                                                     |
| --------------- | ------------------------------------------------------------ |
| `match()`       | 确定正则是否从字符串的开头匹配。                             |
| `search()`      | 扫描字符串，查找此正则匹配的任何位置。                       |
| `findall()`     | 找到正则匹配的所有子字符串，并将它们作为列表返回。           |
| `finditer()`    | 找到正则匹配的所有子字符串，并将它们返回为一个 [iterator](https://docs.python.org/zh-cn/3/glossary.html#term-iterator)。 |

## 常用匹配字符

**元字符**

`. ^ $ * + ? { } [ ] \ | ( )`

1. `.` **(点号)**

- **含义**: 匹配除换行符(`\n`)之外的任何单个字符。
- **示例**:

```python
import re
print(re.findall(r'a.b', 'a1b a2b acb a\nb'))  # 输出 ['a1b', 'a2b', 'acb']
```

2. `^` **(脱字符)**

- **含义**: 匹配字符串的开始位置。
- **示例**:

```python
print(re.findall(r'^Hello', 'Hello World'))  # 输出 ['Hello']
print(re.findall(r'^Hello', 'Hi Hello World'))  # 输出 []
```

3. `$` **(美元符)**

- **含义**: 匹配字符串的结束位置。
- **示例**:

```python
print(re.findall(r'World$', 'Hello World'))  # 输出 ['World']
print(re.findall(r'World$', 'World Hello'))  # 输出 []
```

4. `*` **(星号)**

- **含义**: 匹配前面的字符零次或多次。
- **示例**:

```python
print(re.findall(r'ab*', 'a ab abb abbb'))  # 输出 ['a', 'ab', 'abb', 'abbb']
```

5. `+` **(加号)**

- **含义**: 匹配前面的字符一次或多次。
- **示例**:

```python
print(re.findall(r'ab+', 'a ab abb abbb'))  # 输出 ['ab', 'abb', 'abbb']
```

6. `?` **(问号)**

- **含义**: 匹配前面的字符零次或一次（非贪婪匹配）。
- **示例**:

```python
print(re.findall(r'ab?', 'a ab abb abbb'))  # 输出 ['a', 'ab', 'ab', 'ab']
```

7. `{}` **(大括号)**

- **含义**: 匹配前面的字符指定次数。
- **语法**: `{m}`表示匹配前一个字符`m`次，`{m,n}`表示匹配前一个字符至少`m`次，至多`n`次。
- **示例**:

```python
print(re.findall(r'ab{2}', 'a ab abb abbb'))  # 输出 ['abb']
print(re.findall(r'ab{1,3}', 'a ab abb abbb'))  # 输出 ['ab', 'abb', 'abbb']
```

8. `[]` **(中括号)**

- **含义**: 匹配方括号内的任意一个字符，`-`可以用于指定字符范围。
- **示例**:

```python
print(re.findall(r'[aeiou]', 'hello world'))  # 输出 ['e', 'o', 'o']
print(re.findall(r'[a-z]', 'Hello World 123'))  # 输出 ['e', 'l', 'l', 'o', 'o', 'r', 'l', 'd']
```

9. `|` **(竖线)**

- **含义**: 表示逻辑或，匹配竖线两边的任意一个模式。
- **示例**:

```python
print(re.findall(r'cat|dog', 'I have a cat and a dog'))  # 输出 ['cat', 'dog']
```

10. `()` **(小括号)**

- **含义**: 用于分组，可以将部分模式括起来，并将其作为一个整体进行操作。还可以通过分组提取匹配的子串。
- **示例**:

```python
print(re.findall(r'(ab)+', 'ababab ab ab'))  # 输出 ['ab']
```

11. `\` **(反斜杠)**

- **含义**: 用于转义字符。如果一个字符在正则表达式中具有特殊含义，想要匹配这个字符本身，可以在前面加上反斜杠。
- **示例**:

```python
print(re.findall(r'\.', 'a.b.c'))  # 输出 ['.', '.']
```

12. `\d`

- **含义**: 匹配任何一个数字，相当于`[0-9]`。
- **示例**:

```python
print(re.findall(r'\d+', 'There are 2 cats and 3 dogs.'))  # 输出 ['2', '3']
```

13. `\D`

- **含义**: 匹配任何一个非数字字符，相当于`[^0-9]`。
- **示例**:

```python
print(re.findall(r'\D+', 'There are 2 cats and 3 dogs.'))  # 输出 ['There are ', ' cats and ', ' dogs.']
```

14. `\w`

- **含义**: 匹配任何一个字母、数字或下划线，相当于`[a-zA-Z0-9_]`。
- **示例**:

```python
print(re.findall(r'\w+', 'This is a test_123.'))  # 输出 ['This', 'is', 'a', 'test_123']
```

15. `\W`

- **含义**: 匹配任何一个非字母、数字或下划线的字符，相当于`[^a-zA-Z0-9_]`。
- **示例**:

```python
print(re.findall(r'\W+', 'This is a test_123.'))  # 输出 [' ', ' ', ' ', '.']
```

16. `\s`

- **含义**: 匹配任何一个空白字符，包括空格、制表符、换页符等，相当于`[ \t\n\r\f\v]`。
- **示例**:

```python
print(re.findall(r'\s+', 'This is a test.'))  # 输出 [' ', ' ', ' ']
```

17. `\S`

- **含义**: 匹配任何一个非空白字符，相当于`[^ \t\n\r\f\v]`。
- **示例**:

```python
print(re.findall(r'\S+', 'This is a test.'))  # 输出 ['This', 'is', 'a', 'test.']
```

18. `\b`

- **含义**: 匹配单词的边界。
- **示例**:

```python
print(re.findall(r'\bword\b', 'word in a sentence'))  # 输出 ['word']
print(re.findall(r'\bword\b', 'password in a sentence'))  # 输出 []
```

19. `\B`

- **含义**: 匹配非单词边界。
- **示例**:

```python
print(re.findall(r'\Bword\B', 'password in a sentence'))  # 输出 ['word']
```

## 常用方法

### 1. 匹配字符串

使用 `re.match()` 函数可以检查一个字符串是否以指定的模式开头。

```python
pattern = r'hello'
text = 'hello world'

match = re.match(pattern, text)
if match:
    print("Found match:", match.group())
else:
    print("No match")
```

- `r'hello'` 是一个原始字符串，用于描述要匹配的模式。
- `re.match()` 尝试从字符串的开头开始匹配模式，如果匹配成功则返回一个匹配对象，否则返回 `None`。

### 2. 搜索字符串

使用 `re.search()` 函数可以在字符串中搜索匹配的子串。

```python
pattern = r'world'
text = 'hello world'

match = re.search(pattern, text)
if match:
    print("Found match:", match.group())
else:
    print("No match")
```

- `re.search()` 在整个字符串中搜索第一个匹配的模式。

### 3. 查找所有匹配项

使用 `re.findall()` 函数可以查找字符串中所有匹配的子串，并返回一个列表。

```python
pattern = r'\d+'
text = 'There are 123 apples and 456 oranges'

matches = re.findall(pattern, text)
print("All matches:", matches)
```

- `r'\d+'` 是一个正则表达式，匹配一个或多个数字。
- `re.findall()` 返回一个包含所有匹配项的列表 `['123', '456']`。

### 4. 替换字符串

使用 `re.sub()` 函数可以替换字符串中匹配的子串。

```python
pattern = r'\d+'
replacement = 'X'
text = 'There are 123 apples and 456 oranges'

new_text = re.sub(pattern, replacement, text)
print("New text:", new_text)
```

- `re.sub()` 将字符串中所有匹配 `r'\d+'` 的子串替换为 `'X'`，输出结果为 `'There are X apples and X oranges'`。

### 5. 分割字符串

使用 `re.split()` 函数可以根据匹配模式分割字符串。

```python
pattern = r'\s+'
text = 'Hello   World'

tokens = re.split(pattern, text)
print("Tokens:", tokens)
```

- `re.split()` 根据模式 `r'\s+'`（一个或多个空白字符）来分割字符串，并返回分割后的列表 `['Hello', 'World']`。

### 注意事项

- 正则表达式中常用的元字符包括 `.`（匹配任意字符）、`*`（零个或多个前面的表达式）、`+`（一个或多个前面的表达式）、`?`（零个或一个前面的表达式）、`^`（匹配字符串的开头）、`$`（匹配字符串的结尾）、`\d`（匹配一个数字字符）、`\w`（匹配字母数字字符或下划线）等。
- 在使用正则表达式时，建议先在在线工具（如 [regex101](https://regex101.com/)）上测试和调试你的正则表达式。

通过这些基本的示例和方法，你可以开始使用 Python 中的 `re` 标准库进行文本匹配、搜索和替换操作。