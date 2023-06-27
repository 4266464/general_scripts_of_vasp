#!/home/phys/qif/anaconda3/bin/python

#import hmac, base64, struct, hashlib, time
import pexpect
import argparse

# def get_hotp_token(secret):
#     k2 = secret.strip().replace(' ','')
#     if len(k2)%8 != 0:
#         k2 += '='*(8-len(k2)%8)
#     key = base64.b32decode(k2, True)
#     msg = struct.pack(">Q", int(time.time())//30)
#     h = bytearray(hmac.new(key, msg, hashlib.sha1).digest())
#     o = h[19] & 15
#     h = str((struct.unpack(">I", h[o:o+4])[0] & 0x7fffffff) % 1000000)
#     return '0'*(6-len(h)) + h if len(h) < 6 else h

def scp(user="ckduan",host="114.214.207.167",password="lIZ@7U$OE%"):
    #scp SC132.vasp ckduan@114.214.207.167:/share/home/ckduan/usr/2sidesniddle/t\
    parser = argparse.ArgumentParser()
    parser.add_argument("files",nargs="+",help="file_source file_target of 'pscp file_source file_target'")
    args = parser.parse_args()
    for i in range(len(args.files)):
        args.files[i]=args.files[i].replace("rr",user+"@"+host+":/share/home/ckduan/usr/2sidesniddle")
    command="scp -r "+" ".join(args.files)
    print(command)
    process = pexpect.spawn(command, timeout=30)
    expect_list = ["password",pexpect.EOF,pexpect.TIMEOUT]
    index = process.expect(expect_list)
    if index == 0:
        process.sendline(password)
        process.read()
    elif index == 1:
        print('EOF')
    else:
        print('TIMEOUT')

scp()
# get_hotp_token("V5WJEI5VLDRFJNROKA2Q4JG2GE")
# scp SC132.vasp ckduan@114.214.207.167:/share/home/ckduan/usr/2sidesniddle/t
# ckduan@114.214.207.167's password:
# lIZ@7U$OE%

# scp p.tar qif@211.86.151.106:/home/phys/qif/2sidesniddle/perovskite/0oh/new/0Cs2NaInCl6/relax_ste_low_pbe0_soc/
# lmz@qif@619b