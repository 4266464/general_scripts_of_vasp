# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 00:59:30 2021

@author: dugue
"""

from __future__ import print_function
import sys
import cv2 as cv

use_mask = False
img = None
templ = None
mask = None
image_window = "Source Image"
result_window = "Result window"
match_method = 0
max_Trackbar = 5

def main(img_path,templ_path,mask_path=False):
    img = cv.imread(img_path, cv.IMREAD_COLOR)
    templ = cv.imread(templ_path, cv.IMREAD_COLOR)
    #mask = cv.imread(mask_path, cv.IMREAD_COLOR)
    
    cv.namedWindow(image_window, cv.WINDOW_AUTOSIZE)
    cv.namedWindow(result_window, cv.WINDOW_AUTOSIZE)
    
    
    #trackbar_label = 'Method: \n 0: SQDIFF \n 1: SQDIFF NORMED \n 2: TM CCORR \n 3: TM CCORR NORMED \n 4: TM COEFF \n 5: TM COEFF NORMED'
    #cv.createTrackbar(trackbar_label, image_window, match_method, max_Trackbar, MatchingMethod)
    
    #MatchingMethod(img,templ,match_method)
    MatchingMethod(img,templ,cv.TM_CCORR_NORMED)
    
    #cv.waitKey(0)
    return 0
    
def MatchingMethod(img,templ,param):
    global match_method
    match_method = param
    
    img_display = img.copy()
    
    #method_accepts_mask = (cv.TM_SQDIFF == match_method or match_method == cv.TM_CCORR_NORMED)
    #result = cv.matchTemplate(img, templ, match_method, None, mask)
    result = cv.matchTemplate(img, templ, match_method)
    cv.normalize(result, result, 0, 1, cv.NORM_MINMAX, -1)
    _minVal, _maxVal, minLoc, maxLoc = cv.minMaxLoc(result, None)
    
    print(_minVal, _maxVal, minLoc, maxLoc)
    # if (match_method == cv.TM_SQDIFF or match_method == cv.TM_SQDIFF_NORMED):
    #     matchLoc = minLoc
    # else:
    #     matchLoc = maxLoc
    
    
    #cv.rectangle(img_display, matchLoc, (matchLoc[0] + templ.shape[0], matchLoc[1] + templ.shape[1]), (0,0,0), 2, 8, 0)
    #cv.rectangle(result, matchLoc, (matchLoc[0] + templ.shape[0], matchLoc[1] + templ.shape[1]), (0,0,0), 2, 8, 0)
    #cv.imshow(image_window, img_display)
    #cv.imshow(result_window, result)
    
    pass

main("1628957939581.jpg","tpl1628823344588.png")

















# use_mask = False
# img = None
# templ = None
# mask = None
# image_window = "Source Image"
# result_window = "Result window"
# match_method = 0
# max_Trackbar = 5

# def main(img_path,templ_path,mask_path=False):
#     img = cv.imread(img_path, cv.IMREAD_COLOR)
#     templ = cv.imread(templ_path, cv.IMREAD_COLOR)
#     #mask = cv.imread(mask_path, cv.IMREAD_COLOR)
    
#     cv.namedWindow(image_window, cv.WINDOW_AUTOSIZE)
#     cv.namedWindow(result_window, cv.WINDOW_AUTOSIZE)
    
    
#     trackbar_label = 'Method: \n 0: SQDIFF \n 1: SQDIFF NORMED \n 2: TM CCORR \n 3: TM CCORR NORMED \n 4: TM COEFF \n 5: TM COEFF NORMED'
#     cv.createTrackbar(trackbar_label, image_window, match_method, max_Trackbar, MatchingMethod)
    
#     MatchingMethod(img,templ,match_method)
    
#     cv.waitKey(0)
#     return 0
    
# def MatchingMethod(img,templ,param):
#     global match_method
#     match_method = param
    
#     img_display = img.copy()
    
#     method_accepts_mask = (cv.TM_SQDIFF == match_method or match_method == cv.TM_CCORR_NORMED)
#     if (use_mask and method_accepts_mask):
#         result = cv.matchTemplate(img, templ, match_method, None, mask)
#     else:
#         result = cv.matchTemplate(img, templ, match_method)
    
    
#     cv.normalize(result, result, 0, 1, cv.NORM_MINMAX, -1)
    
#     _minVal, _maxVal, minLoc, maxLoc = cv.minMaxLoc(result, None)
    
    
#     if (match_method == cv.TM_SQDIFF or match_method == cv.TM_SQDIFF_NORMED):
#         matchLoc = minLoc
#     else:
#         matchLoc = maxLoc
    
    
#     cv.rectangle(img_display, matchLoc, (matchLoc[0] + templ.shape[0], matchLoc[1] + templ.shape[1]), (0,0,0), 2, 8, 0)
#     cv.rectangle(result, matchLoc, (matchLoc[0] + templ.shape[0], matchLoc[1] + templ.shape[1]), (0,0,0), 2, 8, 0)
#     cv.imshow(image_window, img_display)
#     cv.imshow(result_window, result)
    
#     pass

# main("pot_1Cs2HfCl6.png","pot_1Cs2HfCl6.png")