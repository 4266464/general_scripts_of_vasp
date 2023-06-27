# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 02:10:20 2021

@author: dugue
"""

import win32api, win32con, win32gui

# hwnd_title={}
# def get_all_hwnd(hwnd,mouse):
#     if win32gui.IsWindow(hwnd) and win32gui.IsWindowEnabled(hwnd) and win32gui.IsWindowVisible(hwnd):
#         hwnd_title[hwnd]=win32gui.GetWindowText(hwnd)
# win32gui.EnumWindows(get_all_hwnd,0)
 
# for h,t in hwnd_title.items():
#     if not t=="":
#         print(h, t)
        
# def get_window_pos(name):
#     name = name
#     handle = win32gui.FindWindow(0, name)
#     if handle == 0:
#         return None
#     else:
#         return win32gui.GetWindowRect(handle)
# x1, y1, x2, y2 = get_window_pos("《战舰世界》")
# print(x1,y1,x2,y2)

import time
import numpy as np
from PIL import ImageGrab
# 每抓取一次屏幕需要的时间约为1s,如果图像尺寸小一些效率就会高一些
img = ImageGrab.grab(bbox=(-1800, 213, -159, 972))
img.show()
#img = np.array(img.getdata(), np.uint8).reshape(img.size[1], img.size[0], 3)


# PIL.Image转换成OpenCV格式：

#     import cv2  
#     from PIL import Image  
#     import numpy  
      
#     image = Image.open("plane.jpg")  
#     image.show()  
#     img = cv2.cvtColor(numpy.asarray(image),cv2.COLOR_RGB2BGR)  
#     cv2.imshow("OpenCV",img)  
#     cv2.waitKey()  
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# cv2保存图片用cv2.imwrite("/home/1.jpg" ,frame * 1) # *1才为彩色,后面*1可以不写
# cv2看图片大小用 img.shape # 它的输出是(480, 640, 3)，记住这里宽是480,长是640，深度是3色的彩色
# cv2裁剪用img = img[60:420, 60:580, :] #eg:[宽为60~（480-60），长为60～（640-60),第三个是选择全部深度]

# OpenCV转换成PIL.Image格式：

#     import cv2  
#     from PIL import Image  
#     import numpy  
      
#     img = cv2.imread("plane.jpg")  
#     cv2.imshow("OpenCV",img)  
#     image = Image.fromarray(cv2.cvtColor(img,cv2.COLOR_BGR2RGB))  
#     image.show()  
#     cv2.waitKey() 
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# PIL的Image保存图片用img.save(“001.jpg”)
# PIL的Image看图片大小用frame.size
# PIL的Image裁剪用crop

# # box = (100, 100,lab_w-100,lab_h-100)  # 左 上 右 下（控制裁剪框大小）(前两个小于后两个数值)
# # img = img.crop(box)
# 1
# 2
# 还有tabel里面画框用cv2.rectangle

# # 画框-->cv2.rectangle
# frame = cv2.rectangle(frame, (100,100 ), (lab_w-100,lab_h-100 ), (0, 0, 255), 6)  #后面两个分别是颜色和框的宽度
# 1
# 2
# 不要cv2转换PIL的Image，此时还是480宽，640高，再来放大原来label大小再来crop裁剪，这样放大会丢失所以还是原图修改，直接按比例就好附上我自己的代码来提高图像稳定性

#     def Preheat(self):
#         temperature = True  # 用来对比，有按下后没按下两者稳定性是否有很大差别，True表示预热去掉边框外的影响，False表示原来整张图加载进去的
#         capture = self.capture  # 当然我前面def __init__里面有self.capture = cv2.VideoCapture(0)
#         # 获取一帧
#         ret,frame = capture.read()
#         # print(frame.shape)
#         # print("lab_w:%d" % (lab_w))
#         # print("lab_h:%d" % (lab_h))

#         frame = cv2.flip(frame, flipCode=1)  # 水平翻转，我label是1356长,837宽
#         x = int(100/837*480)  # 因为下面不准有小数所以转int，有那么一行误差啦不影响
#         w = int(480-100/837*480)
#         y = int(100/1356*640)
#         h = int(640-100/1356*640)
#         frame = frame[x:w, y:h, :]  # 这个cv2里面的裁剪，eg:[宽为60~（480-60），长为60～（640-60),第三个是选择全部深度]
#         cv2.imwrite("000.jpg",frame)
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# 10
# 11
# 12
# 13
# 14
# 15
# 16
# 还有提醒一下：如果是显示视频的话，还是用cv2.imshow(‘frame’,img)，因为死循环while image.show()会不断创建新进程。


# 1.获取窗口左上角及右下角坐标
# import win32api, win32con, win32gui
# def get_window_pos(name):
#     name = name
#     handle = win32gui.FindWindow(0, name)
#  # 获取窗口句柄
#  if handle == 0:
#  return None
#  else:
#  return win32gui.GetWindowRect(handle)
# x1, y1, x2, y2 = get_window_pos('暴雪战网')
# print(x1,y1,x2,y2)
# 结果：

# F:\push\20190929>python 1.py
# (349, 83, 1549, 1013)
# 其中窗口信息(x1, y1, x2, y2)，(x1, y1)是窗口左上角的坐标，(x2, y2)是窗口右下角的坐标。我们可以利用这个信息配合PIL进行截图。但是在这之前，我们还要解决两个问题：

# 该窗口并不在当前的界面上，被其他的软件覆盖到底层中，这时候需要高亮窗口。
# 该窗口被最小化怎么办？
# 2.win32gui 高亮窗口
# 为了使得被叠在底层的窗口能放到最上层显示，我们需要拿到窗口的handle，对其执行高亮操作，其实很简单，我们刚刚获得坐标信息的时候已经得到handle了，只需要做一下简单的更改即可。

# import win32api, win32con, win32gui
# def get_window_pos(name):
#     name = name
#     handle = win32gui.FindWindow(0, name)
#  # 获取窗口句柄
#  if handle == 0:
#  return None
#  else:
#  # 返回坐标值和handle
#  return win32gui.GetWindowRect(handle), handle
# (x1, y1, x2, y2), handle = get_window_pos('暴雪战网')
# text = win32gui.SetForegroundWindow(handle)
# 这样就能将被覆盖到底层的窗口放到最上层，如下图所示。


# 3. 还原最小化窗口
# 还有一种特殊情况就是窗口被缩小了，这时候我们就需要还原最小化窗口，其实也非常简单，只要利用win32gui和win32con向该窗口发送一个信息即可。

# import win32api, win32con, win32gui
# def get_window_pos(name):
#     name = name
#     handle = win32gui.FindWindow(0, name)
#  # 获取窗口句柄
#  if handle == 0:
#  return None
#  else:
#  # 返回坐标值和handle
#  return win32gui.GetWindowRect(handle), handle
# (x1, y1, x2, y2), handle = get_window_pos('暴雪战网')
# win32gui.SendMessage(handle, win32con.WM_SYSCOMMAND, win32con.SC_RESTORE, 0)
# # 发送还原最小化窗口的信息
# win32gui.SetForegroundWindow(handle)
# # 设为高亮 
# 效果如图所示：


# 4.截图
# 有了PIL模块和窗口的坐标后，我们想截图可非常简单。PIL 模块安装：
# pip install pillow
# 安装完就可以试一下我们的完整代码了，如下：

# import win32api, win32con, win32gui
 
# def get_window_pos(name):
#     name = name
#     handle = win32gui.FindWindow(0, name)
#  # 获取窗口句柄
#  if handle == 0:
#  return None
#  else:
#  # 返回坐标值和handle
#  return win32gui.GetWindowRect(handle), handle
# (x1, y1, x2, y2), handle = get_window_pos('暴雪战网')
# win32gui.SendMessage(handle, win32con.WM_SYSCOMMAND, win32con.SC_RESTORE, 0)
# # 发送还原最小化窗口的信息
# win32gui.SetForegroundWindow(handle)
# # 设为高亮
# from PIL import Image, ImageGrab
# img_ready = ImageGrab.grab((x1, y1, x2, y2))
# # 截图
# img_ready.show()
# # 展示 