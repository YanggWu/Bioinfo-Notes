# Git

## 配置相关账户

```bash
git config --global user.name "Your Name"
git config --global user.email "email@example.com"
```

## 基本使用

创建一个新的本地仓库

```bash
# 1. 初始化 Git 仓库

git init	# 在当前目录下生成一个 .git 文件夹，用于存储版本历史、配置信息和对象等。
```

## 远程仓库相关

```bash
# 1.查看当前远程仓库
git remote -v
# 2.移除原始远程仓库（如果存在不相关的仓库）
git remote remove origin

# 3.添加自己的 GitHub 仓库作为新的远程仓库
git remote add origin https://github.com/your-username/nvim.git	# 使用 HTTPS
git remote add origin git@github.com:your-username/nvim.git		# 使用 SSH
git remote -v	# 确认新仓库已正确关联

# 4.更新remote地址
git remote set-url origin xxx
# 例如：修改用户名后需要更新地址
git remote set-url origin https://github.com/yanggwu/Record.git

# 5.推送更改到远程仓库

```

## 常用命令

```bash
git init
git status
git add
git commit
git commit -m "Update"
touch .gitignore  #忽略的文件
git branch new-branch # 创建新分支
git branch # 查看当前分支
git checkout -b new-branch # 创建一个新的分支
git checkout new-branch #切换分支
git commit -a -m "Update" #快速从工作区提交到本地仓库
git branch -d new-branch # 或者 -D 直接删除分支
git merge new-branch # 将该分支合并到当前分支
git clone hppts://xxxxx
git pull
```



## 常见问题

**Failed to connect to github.com port 443 解决方案**

1. 尝试取消代理

```shell
git config --global --unset http.proxy
git config --global --unset https.proxy
```

2. 配置http代理 
有socks5和http两种协议，有使用的代理软件决定，不确定可以都试一下。

主机号和端口号可以从代理软件上查找。q

```bash
# 127.0.0.1 是一个特殊的 IP 地址，指的是本地计算机，也称为 localhost。
git config --global http.proxy 127.0.0.1:7890
git config --global https.proxy 127.0.0.1:7890

git config --global http.proxy socks5 127.0.0.1:7890
git config --global https.proxy socks5 127.0.0.1:7890
```

**push 冲突**
推送冲突的原因一般是远程仓库发生改变，提交者的版本库小于远程仓库。

可以先 git pull 实现同步之后，在 git push 推送到远程仓库。
