# Pandas

## 文件读取

Pandas 支持多种格式的数据文件，如 CSV、Excel、JSON、HTML 和 SQL 数据库等。最常用的读取数据的函数是 `read_csv()` 和 `read_excel()`。

pd.read_csv

```py
pd.read_csv(filepath_or_buffer, 
			sep=',', header='infer', 
			names=None, index_col=None, 
			usecols=None, dtype=None, 
			nrows=None, skiprows=None, 
			na_values=None, encoding=None)

```

`sep` 指定分隔符，默认为制表符“ ,”

`header` 参数控列名所在的行。默认值：`"infer"`自动推断文件的头部行号。

`names` 参数用于手动指定列名，尤其是在文件没有列名时

`index_col` 指定哪些列应被用作行索引。

`usecols` 控制从文件中读取哪些列。

`dtype` 指定列的数据类型,可以使用 dict 为每列指定不同的数据类型。

```py
# 使用逗号作为分隔符 (CSV 文件)
df = pd.read_csv('file.csv', sep=',')

# 使用制表符作为分隔符 (TSV 文件)
df = pd.read_csv('file.tsv', sep='\t')
```

## 数据索引

提供了多种灵活的方式来索引和选取数据，以下是几种常用的索引方式

=== "df.column"
	通过点运算符直接访问 DataFrame 的列，类似于访问对象的属性。适用于列名为合法的Python变量时，即不包含空格、特殊字符等	
	```py
	df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
	df.A  # 相当于 df['A']
	```	
=== "df[ ]"
	df[ ] 使用 `[]` 来选取 DataFrame 中的一列或多列，也可以通过切片和条件表达式过滤行。
	
	```py
	df[:4] # 切片的方式选取前四行
	df[df.value > 80] # 布尔型数组（过滤行）
	df["education"] # 单个值选取列
	df[["education", "salary"]] # 序列选取多列
	```

=== "df.loc"
	基于标签（`label`）的索引方式，用来索引行、列。可以通过标签同时索引行和列，灵活性较高。
	
	```py
	# 按行索引, 注意df的index是数字标签，还是字符串标签。使用的标签一定是df中实际存在的
	df.loc[0]  		# 选取index为 0 的行
	df.loc[0, 'A']  # 选取index为 0 的行， 'A' 列的值
	
	```

=== "df.iloc"
	与df.loc 类似，但是基于位置（`position`）的索引方式，通过整数索引行和列。
	```py
	df.iloc[0]  # 选取第 0 行
	df.iloc[0, 1]  # 选取第 0 行第 1 列的值
	df.iloc[0:2, 0:2]  # 选取第 0 到第 2 行，第 0 到第 2 列
	```
	注意理解df.loc 和df.iloc的区别

## 数据过滤

布尔条件过滤 

使用布尔条件结合数据索引工具过滤 DataFrame 中的行

```py
df[df['A'] > 2]  # 选取 'A' 列大于 2 的所有行
iris[iris.species == 'setosa']
```

**多条件过滤**： 使用逻辑运算符（`&` 和 `|`）进行多个条件的组合。**注意**：在多个条件时，需用括号将每个条件括起来。

```py
# 选取 'A' 大于 2 且 'B' 小于 5 的行
df[(df['A'] > 2) & (df['B'] < 5)]

# 选取 'A' 小于 2 或 'B' 大于 5 的行
df[(df['A'] < 2) | (df['B'] > 5)]
```



### `.query()`

使用查询语法过滤 DataFrame，类似 SQL 语句。

```
df.query('A > 2 and B < 5')  # 选取 'A'列 大于 2 且 'B'列 小于 5 的行
```



## 常用函数

### apply函数

用于在 DataFrame 或 Series 上应用自定义函数或内置函数，以对数据进行元素级或向量级的操作。

- **参数**：
  - func：要应用于每个元素的函数。可以是自定义函数、内置函数或 lambda 函数。
  - axis：指定应用函数的轴方向。默认为0，表示沿着每列应用函数；1 表示沿着每行应用函数。
  - args 和 kwargs：可选参数，用于传递给函数的额外参数和关键字参数。