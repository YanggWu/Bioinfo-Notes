# SSH登陆服务器

SSH登录到远程服务器，可以有两种方法，一种是使用户名和用户密码，另外一种方式是SSH密钥，因为使用用户密码的登录方式相对密钥登录的方式容易被暴力破解，因为推荐使用密钥登录的方式登录到远程服务器。

## Terminus

推荐使用[Terminus](https://www.termius.com/) 登陆远程服务器，支持全平台和数据传输。

## vscode

在config文件中添加 KbdInteractiveAuthentication yes 以使用动态口令登陆。

```txt
Host cluster
    HostName login_ip
    User username
    Port login_port
    KbdInteractiveAuthentication yes
```
