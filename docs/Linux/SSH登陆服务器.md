# SSH登陆服务器

SSH 登录到远程服务器，可以有两种方法，一种是使用户名和用户密码，另外一种方式是SSH密钥，因为使用用户密码的登录方式相对密钥登录的方式容易被暴力破解，因为推荐使用密钥登录的方式登录到远程服务器。

如果是个人本地服务器，通常是ubuntu系统。需要先配置SSH

```bash
sudo apt update
sudo apt install openssh-server	# 安装完成后，SSH服务会自动启动

# 检查SSH服务器的状态
sudo systemctl status ssh

# 确定机器 iP 地址。
ip a
```

## Xshell
Xshell 是一款 Windows 下的流行老牌 SSH 客户端，正版软件为商业收费软件，一般用户可以使用免费版本。

<div class="grid cards" markdown>

- :simple-shelly: XSHELL 官网 [:octicons-arrow-right-24: <a href="https://www.netsarang.com/en/xshell/" target="_blank"> 传送门 </a>](#)

</div>
    

## Terminus

Terminus 是一款跨平台的终端工具，支持 SSH、Telnet 等协议，具备简洁美观的界面，并且支持插件扩展。

<div class="grid cards" markdown>

- :octicons-terminal-16: Terminus 官网 [:octicons-arrow-right-24: <a href="https://www.termius.com/" target="_blank"> 传送门 </a>](#)

</div>

!!! note
    个人推荐 terminus 登陆远程服务器，支持全平台和数据传输。

## VSCode

使用 Visual Studio Code（VSCode） 连接远程服务器可以方便地进行远程开发，尤其是在生物信息学或其他领域中，需要处理远程高性能计算集群时，VSCode 提供了友好的开发环境。通过 VSCode 用户可以在本地编辑、运行和调试远程服务器上的代码。

首先安装 `Remote - SSH`, `Remote Development` 等相关插件

=== "密码登陆"
    在远程资源管理器窗口找到设置中的 config 文件进行密码登陆配置。
    ```txt
    Host cluster
        HostName login_ip
        User username
        Port login_port
    ```

    !!! tip
        每次登陆需要在窗口输入密码，推荐使用密钥免密登陆更方便。

=== "密钥登陆"
    首先本地生成密钥，然后将本地生成的 `id_rsa.pub` 公钥的内容复制好。
    远程服务器同样生成密钥。然后创建一个 `authorized_keys` 文件。将本地公钥内容添加进去。
    ```bash
    # 1. 本地生成密钥
    ssh-keygen -t rsa
    # 2. 进入密钥所在目录，复制本地公钥的内容
    cd ~/.ssh
    # 3. 远程服务器同样生成密钥
    ssh-keygen -t rsa
    # 4. 在远程服务器的密钥目录创建 authorized_keys
    vim authorized_keys # 添加本地公钥内容
    ```

=== "动态口令登陆"

    在config文件中添加 KbdInteractiveAuthentication yes 以使用动态口令登陆。
    
    ```txt
    Host cluster
        HostName login_ip
        User username
        Port login_port
        KbdInteractiveAuthentication yes
    ```
