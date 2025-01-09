# Git

Git 是一个分布式版本控制系统，广泛用于软件开发和协作项目中。以下是 Git 的基本操作介绍。

<img src="https://raw.githubusercontent.com/wuayng-owl/Image/main/markdown_image/20231226234933.png" width="650">

## 常用命令

```bash
# 基本使用
git init				  # 在当前目录创建一个新的 Git 仓库
git status			      # 查看工作区和暂存区的状态
git add <文件名>	        # 将更改的文件加入暂存区
git add .		   	      # 添加所有更改的文件
git commit
git commit -m "Update"    # 将暂存区的更改提交到本地仓库
git commit -a -m "Update" # 快速从工作区提交到本地仓库
git pull
touch .gitignore  	      # 忽略的文件

# 分支相关
git branch <分支名>	  	# 创建一个新分支
git branch 				  # 查看当前分支
git checkout <分支名>	    # 切换到指定分支
git checkout -b <分支名> 	# 同时创建并切换到新分支
git branch -d <分支名>  	# 或者 -D 直接删除分支
git merge new-branch 	  # 将该分支合并到当前分支
git clone hppts://xxxxx	  # 从远程仓库复制项目到本地
```

## 配置相关账户

```bash
git config --global user.name "Your Name"
git config --global user.email "email@example.com"
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

### 分支

```bash
# 克隆特定分支
git clone --branch <branch_name> <repository-url>

# 查看远程分支
git branch -r

# 切换到远程分支（本地次无分支时）
git checkout -b feature origin/feature

# 切换到本地已有的分支
git checkout develop

设置本地分支的远程跟踪分支
git branch --set-upstream-to=origin/feature

# 拉取当前分支关联的远程分支更新
git pull
```

## 常见问题

**Failed to connect to github.com port 443 解决方案**

**尝试取消代理**:

```shell
git config --global --unset http.proxy
git config --global --unset https.proxy
```

**配置http代理**: 

有socks5和http两种协议，有使用的代理软件决定，不确定可以都试一下。主机号和端口号可以从代理软件上查找。

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

## 初次运行 Git 前的配置

```bash
# 用户信息配置
git config --global user.name "Wu Yang"
git config --global user.email ywu.info@gmail.com
```

!!! note
    如果用了 `--global` 选项，那么更改的配置文件就是位于你用户主目录下的那个，以后你所有的仓库都会默认使用这里配置的用户信息。如果要在某个特定的仓库中使用其他名字或者电邮，只要去掉 --global 选项重新配置即可，新的设定保存在当前仓库的 .git/config 文件里。

### 配置密钥对

如果是新的服务器上，需要生成 SSH 密钥对

```bash
ssh-keygen -t rsa -b 4096 -C "ywu.info@gmail.com"
```

`-t rsa`：指定使用 RSA 算法生成密钥。

`-b 4096`：设置密钥的长度为 4096 位。

`-C` 后面是用于标识密钥的注释，可以设置为你的邮箱或其他标识信息。

然后会提示你设置密钥文件的保存位置（默认是 `~/.ssh/id_rsa`），以及设置一个 passphrase（密码短语）。你可以选择不设置 passphrase，直接按回车。

**1. 为代码仓库配置SSH密钥:**

生成密钥后，需要将公钥添加到代码仓库 。(`GitHub`, `Gitee`)

```
# 查看公钥：
cat ~/.ssh/id_rsa.pub
```
!!! tips
    登录到你的 GitHub 账户，进入 **Settings** -> **SSH and GPG keys** -> **New SSH key**。

    在 "Title" 中填写一个名称，然后将 `id_rsa.pub` 文件中的内容粘贴到 "Key" 栏中。

**2. 测试 SSH 连接**

```bash
# 1.测试 GitHub
ssh -T git@github.com

# 2.测试 Gitee
ssh -T git@gitee.com
```

### 配置 SSH 客户端

在服务器上配置 SSH 客户端，让它能正确区分使用 GitHub 还是 Gitee 的密钥。

- 编辑 `~/.ssh/config` 文件，如果没有这个文件，可以创建一个。

```bash
vim ~/.ssh/config
```

- 添加以下内容：

```
# GitHub 配置
Host github.com
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_rsa

# Gitee 配置
Host gitee.com
  HostName gitee.com
  User git
  IdentityFile ~/.ssh/id_rsa
```
