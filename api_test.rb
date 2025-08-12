require 'uri'
require 'net/http'
require 'debug'
require 'json'

url = URI.parse("https://jsonplaceholder.typicode.com/posts/1")
response = Net::HTTP.get_response(url)

debugger