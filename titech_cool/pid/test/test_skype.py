#!/usr/bin/env python3

from skpy import Skype
from skpy.msg import SkypeMsg

gmail = 'atlas.jp.itk.pixel.production@gmail.com'

msg = '[titech_cool] this is a test message'
    
password = 'Thus-lcu7-5prA-8og&'
chatid = "19:afd1c48f799c476d95469cbac3ae51ab@thread.skype"

SKPY_DEBUG_HTTP=1
print( 'here' )
sk = Skype(gmail, password) # Tateyama ATLAS Japan ITk                                                                     
print( 'here2' )
sk.chats[chatid].sendMsg(msg, rich=True)

