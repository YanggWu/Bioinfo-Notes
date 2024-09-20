通过搜索Google fonts的可用字体列表，从中可以找到思源宋体。

GOOGLE FONTS 中文字体链接(https://fonts.google.com/?subset=chinese-simplified&script=Hans)

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/image-20240920173533877.png" width="400">

可以获取思源宋体的在线链接`https://fonts.googleapis.com/css2?family=Noto+Serif+SC:wght@200..900&display=swap`

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/202409201737383.png" width="400">

在自定义CSS文件中设置使用该字体

```CSS
/******* 字体、文字相关设置 *******/
body {
    font-family: "Noto Serif SC", system-ui;
    font-weight: 500; /* 设置字体基础字重 */
}
```

 在 mkdocs.yml 文件中配置相关链接。

```yaml
extra_css:
  - stylesheets/main.css
  - stylesheets/extra.css
  - https://fonts.googleapis.com/css2?family=Noto+Serif+SC:wght@200..900&display=swap #自定义思源宋体
```

## 详细自定义字体设置

设置最上方标题大小

```CSS
.md-header__ellipsis .md-header__topic:nth-child(1) .md-ellipsis { /* 最上方标题的样式 */ 
    font-size: 18px;  /* 增大字体大小 */
}
```

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/202409202026499.png" width="400">



设置最上面一排目录样式

```
li.md-tabs__item  {             /* 最上面一级目录的样式 */       
  /* 修改字体 */
  font-weight: 600;
  color: whitesmoke;
}
```

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/202409202030026.png" width="400">