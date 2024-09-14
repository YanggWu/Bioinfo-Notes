# Singularity

**Singularity** 是一种专门用于高性能计算（HPC）和科学计算环境的容器技术，它允许用户在不具有超级用户权限的情况下，创建和运行容器。这使得它非常适合在共享的集群环境中使用，特别是在生物信息学领域，Singularity 是替代 **Docker** 的一个常见工具，尤其是在 HPC 和没有 root 权限的服务器上。

**为什么使用 Singularity？**

1. **无需 root 权限**：Singularity 容器可以在没有超级用户权限的情况下运行，这使其在共享计算资源中尤为重要。
2. **兼容 Docker**：你可以直接使用 Docker 容器镜像，而无需进行任何重大修改。
3. **便于集成 HPC 系统**：许多 HPC 系统支持 Singularity，因为它提供了对文件系统和计算资源的强大支持。
4. **高安全性**：Singularity 容器中用户和主机系统的 UID 是一致的，避免了 Docker 中存在的权限提升问题。

[官方文档](https://www.sylabs.io/guides/3.1/user-guide/)

## 安装

Singularity 需要 **Go** 语言和 **Development Tools**，因此需要先安装这些依赖。

```bash
# 安装Go语言（用于编译 Singularity）
wget https://go.dev/dl/go1.20.5.linux-amd64.tar.gz
sudo tar -C /usr/local -xzf go1.20.5.linux-amd64.tar.gz
echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc
```

从 GitHub 上下载 **Singularity 3.8.7** 的源码

```bash
wget https://github.com/apptainer/singularity/releases/download/v3.8.7/singularity-3.8.7.tar.gz

```



## 使用 Singularity

singularity有许多命令，通过 pull、exec、build 等命令拉取、运行和构建容器。

=== "pull"

    从给定的URL下载容器镜像，常用的有URL有Docker Hub(docker://user/image:tag) 和 Singularity Hub(shub://user/image:tag)
    ```bash
    singularity pull docker://ubuntu:20.04
    ```
    这将会下载一个名为 ubuntu_20.04.sif 的 Singularity 镜像文件

=== "exec"

    运行 Singularity 容器： 使用拉取的镜像，运行 Singularity 容器：
    ```bash
    singularity exec ubuntu_20.04.sif bash
    ```
    这将在容器内部运行 bash，你可以在容器内执行命令。

