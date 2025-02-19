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

### 行筛选

`.query()` 使用查询语法过滤 DataFrame，类似 SQL 语句：

```py
df.query('A > 2 and B < 5')  # 选取 'A'列 大于 2 且 'B'列 小于 5 的行

snp_list = ["vg0100002827", "vg0100002988"]
cultivated_df.query('SNP in @snp_list')		# 使用列表内的元素进行过滤
```

## 列相关操作

**重命名列**

可以通过 `rename()` 函数只修改特定的列名：

```py
# 创建一个示例 DataFrame
data = {'A': [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]}
df = pd.DataFrame(data)

# 使用 rename 修改特定列名
df.rename(columns={'A': 'Alpha', 'B': 'Beta'}, inplace=True)

print(df)
```

- `columns={'A': 'Alpha', 'B': 'Beta'}` 指定了要修改的列及其对应的新名称。

- `inplace=True` 会直接在原 `DataFrame` 上修改。

如果你只知道列的位置（索引）而不想列出所有列名，可以使用列索引来修改列名：

```py
df.columns.values[0] = 'Alpha'  # 修改第一列
df.columns.values[1] = 'Beta'   # 修改第二列
df.columns.values[2] = 'Gamma'  # 修改第三列

print(df)
```

**筛选列**

可以通过列名筛选你需要的列和调整列顺序：

```py
import seaborn as sns
import pandas as pd
df = sns.load_dataset('iris')

# 选择列
df_filtered = df[['sepal_length', 'species']]

# 调整列顺序，只需重新指定列顺序即可。
df_filtered = df_filtered[['species', 'sepal_length']]

# 使用 loc 和 iloc 也可以选择列
df = df.loc[:, ['sepal_length', 'species']]
df = df.iloc[:, [0, 4]]
```

根据列索引排序调整列顺序:

```py
df.sort_index(axis=1, ascending=False)
```

根据条件来筛选列， 比如，选择数据类型为数值类型的列：

```py
# 筛选出数值类型的列
df_numeric = df.select_dtypes(include=['number'])

# 筛选出字符型列
df_string = df.select_dtypes(include=['object'])
```

根据列名包含某些特定字符进行筛选:

```py
# 筛选列名中包含"gene"的列
df_filtered = df.loc[:, df.columns.str.contains('gene')]
```

## 常用函数

### apply函数

用于对 DataFrame 或 Series 进行逐元素操作，或应用一个函数来进行行或列级别的计算。它能够将用户定义的函数应用到 DataFrame 或 Series 的每一行或每一列，或者每个元素，灵活性非常高。

- func：要应用于每个元素的函数。可以是自定义函数、内置函数或 lambda 函数。
- axis：指定应用函数的轴方向。默认为0，表示沿着每列应用函数；1 表示沿着每行应用函数。
- args 和 kwargs：可选参数，用于传递给函数的额外参数和关键字参数。

Series 的 `apply` 函数:

```py
# 创建 Series
s = pd.Series([1, 2, 3, 4, 5])

# 使用 apply 对每个元素进行平方操作
squared = s.apply(lambda x: x ** 2)

# 应用自定义函数
def custom_function(x):
	return x + 10
result = s.apply(custom_function)
```

