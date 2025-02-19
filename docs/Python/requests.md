# requests

`requests` 是一个用于简化 HTTP 请求的第三方 Python 库，广泛用于发送 HTTP 请求（如 GET、POST、PUT、DELETE 等），并处理与服务器进行通信的过程。它使得向 Web 服务器发送请求和处理响应数据变得非常简单。

**主要功能：**

**发送 HTTP 请求**：能够发送各种类型的 HTTP 请求，如 `GET`、`POST`、`PUT`、`DELETE` 等。

**处理响应**：能够处理来自服务器的响应，包括响应内容、状态码、响应头等。

**自动处理 Cookie**：自动管理和发送 Cookie。

**处理 JSON**：轻松处理 JSON 数据。

**支持会话**：可以通过会话对象保持与服务器的持续连接。

**文件上传和下载**：支持文件的上传和下载。

**支持 SSL**：可以通过 HTTPS 协议发送请求。

**异常处理**：提供丰富的异常处理机制，方便捕获网络请求中的错误。

### 常用功能示例：

1. **发送 GET 请求**：

   ```py
   import requests
   
   response = requests.get('https://api.github.com')
   print(response.status_code)  # 返回状态码
   print(response.json())  # 如果返回是JSON格式，可以直接调用 json() 方法解析
   ```

2. **发送 POST 请求**：

   ```py
   data = {'key': 'value'}
   response = requests.post('https://httpbin.org/post', data=data)
   print(response.json())  # 打印返回的 JSON 数据
   ```

3. **处理 JSON 数据**：

   ```py
   response = requests.get('https://jsonplaceholder.typicode.com/posts')
   posts = response.json()  # 解析 JSON 数据
   print(posts[0]['title'])
   ```

4. **上传文件**：

   ```py
   files = {'file': open('file.txt', 'rb')}
   response = requests.post('https://httpbin.org/post', files=files)
   print(response.json())
   ```

5. **使用会话**： 会话对象可以在多个请求之间保持连接，从而避免重复的连接和认证过程。

   ```py
   with requests.Session() as session:
       session.auth = ('user', 'pass')
       response = session.get('https://httpbin.org/basic-auth/user/pass')
       print(response.status_code)
   ```

6. **处理异常**： `requests` 提供了异常机制，可以帮助我们处理网络请求中的错误。

   ```py
   try:
       response = requests.get('https://nonexistent.url')
       response.raise_for_status()  # 如果返回的状态码是 4xx 或 5xx，会抛出 HTTPError 异常
   except requests.exceptions.RequestException as e:
       print(f"请求错误: {e}")
   ```
