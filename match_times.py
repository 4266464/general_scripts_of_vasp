import re
# 不定分组匹配
def match_times(sep, char, string):
    # char="\w" for element "\d" for nelement
    # 保持格式一致
    if sep == "\s":
        string = " " + string
    else:
        string = sep + string
    match_result = re.search("(" + sep + "+(" + char + "+))", string)
    match_list = []
    while match_result:
        unit = match_result.group(2)
        if char == "\d":
            unit = int(unit)
        match_list.append(unit)
        match_result = re.search("(" + sep + "+(" + char + "+)){" + str(len(match_list) + 1) + "}", string)
    return match_list