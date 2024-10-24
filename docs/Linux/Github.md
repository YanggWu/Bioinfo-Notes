# Git

## 常见报错的解决方法

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
