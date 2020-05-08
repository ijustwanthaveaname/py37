from urllib.request import urlopen
import json
json_url='http://raw.githubusercontent.com/muxuezi/btc/master/btc_close_2017.json'
response=urlopen(json_url)
#读取数据
req=response.read()
#将数据写入文件
with open('btc_close_2017_urllib.json','wb') as f:
    f.write(req)
#加载json格式
file_urllib=json.loads(req)
print(file_urllib)

#或者用requests模块下载，更简单
import requests
json_url='http://raw.githubusercontent.com/muxuezi/btc/master/btc_close_2017.json'
req=requests.get(json_url)
#将数据写入文件
with open('btc_close_2017_request.json','w') as f:
    f.write(req.text)
file_requests=req.json()
