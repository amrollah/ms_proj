import httplib2
h = httplib2.Http(".cache")
print('gggg')
(resp_headers, content) = h.request("http://google.com/", "GET")
print(resp_headers)
print(content)