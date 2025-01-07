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

## Gitee 使用示例

### 初次运行 Git 前的配置

```bash
# 用户信息配置
git config --global user.name "Wu Yang"
git config --global user.email ywu.info@gmail.com
```

!!! note
    如果用了 `--global` 选项，那么更改的配置文件就是位于你用户主目录下的那个，以后你所有的仓库都会默认使用这里配置的用户信息。如果要在某个特定的仓库中使用其他名字或者电邮，只要去掉 --global 选项重新配置即可，新的设定保存在当前仓库的 .git/config 文件里。
