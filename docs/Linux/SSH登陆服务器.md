# SSH登陆服务器

SSH登录到远程服务器，可以有两种方法，一种是使用户名和用户密码，另外一种方式是SSH密钥，因为使用用户密码的登录方式相对密钥登录的方式容易被暴力破解，因为推荐使用密钥登录的方式登录到远程服务器。

## Terminus

推荐使用[Terminus](https://www.termius.com/) 登陆远程服务器，支持全平台和数据传输。

## vscode

使用 Visual Studio Code（VSCode） 连接远程服务器可以方便地进行远程开发，尤其是在生物信息学或其他领域中，需要处理远程高性能计算集群时，VSCode 提供了友好的开发环境。通过 VSCode 的 Remote Development 插件，用户可以在本地编辑、运行和调试远程服务器上的代码。

### 密码登陆

### 密钥登陆

### 动态口令登陆

在config文件中添加 KbdInteractiveAuthentication yes 以使用动态口令登陆。

```txt
Host cluster
    HostName login_ip
    User username
    Port login_port
    KbdInteractiveAuthentication yes
```
