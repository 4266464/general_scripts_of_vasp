def flatten(items):
    for x in items:
        # 终止条件，检验是否为可迭代对象
        if hasattr(x,'__iter__') and not isinstance(x, (str, bytes)):
            #Python2写法
            # for sub_x in flatten(x):
            # yield sub_x
            #Python3写法
            yield from flatten(x)
        else:
            yield x
le = list(flatten(list))