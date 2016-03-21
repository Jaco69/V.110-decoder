#!/usr/bin/env python
# https://github.com/Jaco69/V.110-decoder
from sys import argv
from sys import exit
from os.path import splitext
import datetime
import string

def read_byte(f):
  while True:
    b = f.read(1)
    if b:
      yield ord(b)
    else:
      break

def read_t_UI8(f):
  return next(read_byte(f))

def read_t_UI16(f):
  return next(read_byte(f)) + (next(read_byte(f)) << 8)

def read_t_UI32(f):
  return next(read_byte(f)) + (next(read_byte(f)) << 8) + (next(read_byte(f)) << 16) + (next(read_byte(f)) << 24)

def read_t_UI64(f):
  return (next(read_byte(f)) + (next(read_byte(f)) << 8) + (next(read_byte(f)) << 16) + (next(read_byte(f)) << 24)
          + (next(read_byte(f)) << 32) + (next(read_byte(f)) << 40) + (next(read_byte(f)) << 48) + (next(read_byte(f)) << 56))


def hdlc(hdlcbytes, n):
  if bytestream_max_end[n] < bytestream_end[n] + len(hdlcbytes):
    for bytestream_max_end[n] in range(bytestream_max_end[n], bytestream_end[n] + len(hdlcbytes)):
      bytestream[n].append(0)
  for i in range(0, len(v110bytes)):
    bytestream[n][bytestream_end[n]] = hdlcbytes[i]  
    stream_end[n] += 1
  # hdlc frames start and end with 7E
  # so throw away everything before the first 7E
  while bytestream[n][0] != 0x7E:
    for i in range(0, bytestream_end[n]):
      bytestream[n][i] = bytestream[n][i]+1
    bytestream_end[n] -= 1
  # check for second 7E
  for i in range(1, bytestream_end[n]+1):
    if bytestream[n][i] == 0x7E:
      break;
      
def check_bitstream(b, n):
  if (bit_rate%4800 == 0):
    if bitstream_max_end[n] < bitstream_end[n] + 48:
      for bitstream_max_end[n] in range(bitstream_max_end[n], bitstream_end[n] + 48):
        bitstream[n].append(0)
    bitstream[n][bitstream_end[n]] = b[1*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[1*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[1*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[1*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[1*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[1*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[2*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[2*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[2*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[2*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[2*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[2*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[3*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[3*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[3*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[3*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[3*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[3*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[4*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[4*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[4*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[4*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[4*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[4*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[6*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[6*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[6*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[6*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[6*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[6*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[7*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[7*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[7*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[7*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[7*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[7*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[8*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[8*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[8*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[8*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[8*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[8*8+7] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[9*8+2] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[9*8+3] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[9*8+4] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[9*8+5] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[9*8+6] ; bitstream_end[n] += 1
    bitstream[n][bitstream_end[n]] = b[9*8+7] ; bitstream_end[n] += 1
    for i in range(bitstream_end[n]-48, bitstream_end[n]):
      out[n].write(str(bitstream[n][i]))
    out[n].write("\t") ; col[n] += 1
  elif (bit_rate == 2400):
    if bitstream_end[n] + 24 > bitstream_max_end[n]:
      for bitstream_max_end[n] in range(bitstream_max_end[n], bitstream_end[n] + 24):
        bitstream[n].append(0)
    bitstream[n][bitstream_end[n]] = b[1*8+2] ; bitstream_end[n] += 1 #b[1*8+3]
    bitstream[n][bitstream_end[n]] = b[1*8+4] ; bitstream_end[n] += 1 #b[1*8+5]
    bitstream[n][bitstream_end[n]] = b[1*8+6] ; bitstream_end[n] += 1 #b[1*8+7]
    bitstream[n][bitstream_end[n]] = b[2*8+2] ; bitstream_end[n] += 1 #b[2*8+3]
    bitstream[n][bitstream_end[n]] = b[2*8+4] ; bitstream_end[n] += 1 #b[2*8+5]
    bitstream[n][bitstream_end[n]] = b[2*8+6] ; bitstream_end[n] += 1 #b[2*8+7]
    bitstream[n][bitstream_end[n]] = b[3*8+2] ; bitstream_end[n] += 1 #b[3*8+3]
    bitstream[n][bitstream_end[n]] = b[3*8+4] ; bitstream_end[n] += 1 #b[3*8+5]
    bitstream[n][bitstream_end[n]] = b[3*8+6] ; bitstream_end[n] += 1 #b[3*8+7]
    bitstream[n][bitstream_end[n]] = b[4*8+2] ; bitstream_end[n] += 1 #b[4*8+3]
    bitstream[n][bitstream_end[n]] = b[4*8+4] ; bitstream_end[n] += 1 #b[4*8+5]
    bitstream[n][bitstream_end[n]] = b[4*8+6] ; bitstream_end[n] += 1 #b[4*8+7]
    bitstream[n][bitstream_end[n]] = b[6*8+2] ; bitstream_end[n] += 1 #b[6*8+3]
    bitstream[n][bitstream_end[n]] = b[6*8+4] ; bitstream_end[n] += 1 #b[6*8+5]
    bitstream[n][bitstream_end[n]] = b[6*8+6] ; bitstream_end[n] += 1 #b[6*8+7]
    bitstream[n][bitstream_end[n]] = b[7*8+2] ; bitstream_end[n] += 1 #b[7*8+3]
    bitstream[n][bitstream_end[n]] = b[7*8+4] ; bitstream_end[n] += 1 #b[7*8+5]
    bitstream[n][bitstream_end[n]] = b[7*8+6] ; bitstream_end[n] += 1 #b[7*8+7]
    bitstream[n][bitstream_end[n]] = b[8*8+2] ; bitstream_end[n] += 1 #b[8*8+3]
    bitstream[n][bitstream_end[n]] = b[8*8+4] ; bitstream_end[n] += 1 #b[8*8+5]
    bitstream[n][bitstream_end[n]] = b[8*8+6] ; bitstream_end[n] += 1 #b[8*8+7]
    bitstream[n][bitstream_end[n]] = b[9*8+2] ; bitstream_end[n] += 1 #b[9*8+3]
    bitstream[n][bitstream_end[n]] = b[9*8+4] ; bitstream_end[n] += 1 #b[9*8+5]
    bitstream[n][bitstream_end[n]] = b[9*8+6] ; bitstream_end[n] += 1 #b[9*8+7]
    for i in range(bitstream_end[n]-24, bitstream_end[n]):
      out[n].write(str(bitstream[n][i]))
    out[n].write("\t") ; col[n] += 1
  elif (bit_rate == 1200):
    if bitstream_end[n] + 12 > bitstream_max_end[n]:
      for bitstream_max_end[n] in range(bitstream_max_end[n], bitstream_end[n] + 12):
        bitstream[n].append(0)
    bitstream[n][bitstream_end[n]] = b[1*8+2] ; bitstream_end[n] += 1 #b[1*8+3] b[1*8+4] b[1*8+5]
    bitstream[n][bitstream_end[n]] = b[1*8+6] ; bitstream_end[n] += 1 #b[1*8+7] b[2*8+2] b[2*8+3]
    bitstream[n][bitstream_end[n]] = b[2*8+4] ; bitstream_end[n] += 1 #b[2*8+5] b[2*8+6] b[2*8+7]
    bitstream[n][bitstream_end[n]] = b[3*8+2] ; bitstream_end[n] += 1 #b[3*8+3] b[3*8+4] b[3*8+5]
    bitstream[n][bitstream_end[n]] = b[3*8+6] ; bitstream_end[n] += 1 #b[3*8+7] b[4*8+2] b[4*8+3]
    bitstream[n][bitstream_end[n]] = b[4*8+4] ; bitstream_end[n] += 1 #b[4*8+5] b[4*8+6] b[4*8+7]
    bitstream[n][bitstream_end[n]] = b[6*8+2] ; bitstream_end[n] += 1 #b[6*8+3] b[6*8+4] b[6*8+5]
    bitstream[n][bitstream_end[n]] = b[6*8+6] ; bitstream_end[n] += 1 #b[6*8+7] b[7*8+2] b[7*8+3]
    bitstream[n][bitstream_end[n]] = b[7*8+4] ; bitstream_end[n] += 1 #b[7*8+5] b[7*8+6] b[7*8+7]
    bitstream[n][bitstream_end[n]] = b[8*8+2] ; bitstream_end[n] += 1 #b[8*8+3] b[8*8+4] b[8*8+5]
    bitstream[n][bitstream_end[n]] = b[8*8+6] ; bitstream_end[n] += 1 #b[8*8+7] b[9*8+2] b[9*8+3]
    bitstream[n][bitstream_end[n]] = b[9*8+4] ; bitstream_end[n] += 1 #b[9*8+5] b[9*8+6] b[9*8+7]
    for i in range(bitstream_end[n]-12, bitstream_end[n]):
      out[n].write(str(bitstream[n][i]))
    out[n].write("\t") ; col[n] += 1
  elif (bit_rate == 600):
    if bitstream_end[n] + 6 > bitstream_max_end[n]:
      for bitstream_max_end[n] in range(bitstream_max_end[n], bitstream_end[n] + 6):
        bitstream[n].append(0)
    bitstream[n][bitstream_end[n]] = b[1*8+2] ; bitstream_end[n] += 1 #b[1*8+3] b[1*8+4] b[1*8+5] b[1*8+6] b[1*8+7] b[2*8+2] b[2*8+3]
    bitstream[n][bitstream_end[n]] = b[2*8+4] ; bitstream_end[n] += 1 #b[2*8+5] b[2*8+6] b[2*8+7] b[3*8+2] b[3*8+3] b[3*8+4] b[3*8+5]
    bitstream[n][bitstream_end[n]] = b[3*8+6] ; bitstream_end[n] += 1 #b[3*8+7] b[4*8+2] b[4*8+3] b[4*8+4] b[4*8+5] b[4*8+6] b[4*8+7]
    bitstream[n][bitstream_end[n]] = b[6*8+2] ; bitstream_end[n] += 1 #b[6*8+3] b[6*8+4] b[6*8+5] b[6*8+6] b[6*8+7] b[7*8+2] b[7*8+3]
    bitstream[n][bitstream_end[n]] = b[7*8+4] ; bitstream_end[n] += 1 #b[7*8+5] b[7*8+6] b[7*8+7] b[8*8+2] b[8*8+3] b[8*8+4] b[8*8+5]
    bitstream[n][bitstream_end[n]] = b[8*8+6] ; bitstream_end[n] += 1 #b[8*8+7] b[9*8+2] b[9*8+3] b[9*8+4] b[9*8+5] b[9*8+6] b[9*8+7]
    for i in range(bitstream_end[n]-6, bitstream_end[n]):
      out[n].write(str(bitstream[n][i]))
    out[n].write("\t") ; col[n] += 1
  # count zeroes, display and store
  # count ones, keep one
  # get byte, display position to stopbit
  if bitstream[n][0] < 0:
    breaks = -bitstream[n][0]
    bit = 1
  else:
    breaks = 0
    bit = 0
  bb = []
  while bit < bitstream_end[n]:
    # count zeroes
    for bit in range(bit, bitstream_end[n]):
      if not bitstream[n][bit]:
        breaks += 1
      else:
        break
    else: # for run out of values 
      bitstream[n][0] = -breaks
      bitstream_end[n] = 1
      break
    if breaks:
      out[n].write(str(breaks) + "bit break\t") ; col[n] += 1
      breaks = 0
    # skip ones  
    for bit in range(bit, bitstream_end[n]):
      if not bitstream[n][bit]:
        break
    else: # for ran out of values
      # keep 1 
      bitstream[n][0] = 1
      bitstream_end[n] = 1
      break
    # we are now at a zero after a one... 
    if bit+9 >= bitstream_end[n]: # first check if there are enough bits to find a byte
      i = 0
      for bit in range(bit-1, bitstream_end[n]): # keep the one
        bitstream[n][i] = bitstream[n][bit] ; i += 1
      bitstream_end[n] = i
      break
    # calculate byte
    databyte = (bitstream[n][bit+1] << 0)\
             + (bitstream[n][bit+2] << 1)\
             + (bitstream[n][bit+3] << 2)\
             + (bitstream[n][bit+4] << 3)\
             + (bitstream[n][bit+5] << 4)\
             + (bitstream[n][bit+6] << 5)\
             + (bitstream[n][bit+7] << 6)\
             + (bitstream[n][bit+8] << 7)
    # is there a stop bit and so a valid byte?
    if bitstream[n][bit+9]:
      # valid so print byte
      bit = bit + 9
      bb.append(databyte)
      out[n].write('{:02X}'.format(databyte) + '\t') ; col[n] += 1
    else:
      # not valid
      if not databyte: # another break
        continue
      else: #invalid byte what to do?
        continue
  if len(bb) > 0:
    for col[n] in range(col[n], 16):
      out[n].write('\t')
    if decode_ascii:  
      out[n].write('"')
      for i in bb:
        if chr(i) in string.printable:
          if chr(i) == '\t':
            out[n].write('\\t')
          elif chr(i) == '\n':  
            out[n].write('\\n')
          elif chr(i) == '\r':  
            out[n].write('\\r')
          elif chr(i) == '\x0b':  
            out[n].write('\\x0b')
          elif chr(i) == '\x0c':  
            out[n].write('\\x0c')
          else:
            out[n].write(chr(i))      
        else:
          out[n].write('\\x' + '{:02x}'.format(i))
      out[n].write('"')
    if decode_hdlc:
      hdlc(bb, n)
  
def check_v110_frame(f, n):
  for col[n] in range(col[n], 3):
    out[n].write('\t')
  col[n] = 3  
  out[n].write('frame\t') ; col[n] += 1
  out[n].write(str(f[1*8+8]) + str(f[3*8+8]) + str(f[6*8+8]) + str(f[8*8+8]))  # S1,S3,S6,S8
  out[n].write('\t') ; col[n] += 1
  out[n].write(str(f[4*8+8]) + str(f[9*8+8]))  # S4,S9
  out[n].write('\t') ; col[n] += 1
  out[n].write(str(f[2*8+8]) + str(f[7*8+8]))  # X2,X7
  out[n].write('\t') ; col[n] += 1
  out[n].write(str(f[5*8+2]) + str(f[5*8+3]) + str(f[5*8+4]) + str(f[5*8+5]) + str(f[5*8+6]) + str(f[5*8+7]) + str(f[5*8+8]))  # 011 E4-E7
  out[n].write('\t') ; col[n] += 1
  check_bitstream(f, n)
  out[n].write('\n') ; col[n] = 1 ; row[n] += 1
    
def adapted_find_v110_frames(bslen, bytestream, n):
  #for ALL baud rates transparant asynchronous
  frame = {}
  #
  # copy new values to pool
  #
  for i in range(0, bslen):
    stream_end[n] += 1
    if stream_end[n] > stream_max_end[n]:
      stream[n].append(bytestream[i])
      stream_max_end[n] += 1 
    else:
      stream[n][stream_end[n]] = bytestream[i]
  # get stored values
  if stream[n][1] < 0:
    nosync = -stream[n][1]
    framestart = 2
  else:
    nosync = 0
    framestart = 1
  while framestart <= stream_end[n]:
    # check stream for flush
    if (stream[n][framestart] == flush):
      if nosync:
        for col[n] in range(col[n], 3):
          out[n].write('\t')
        out[n].write(str(nosync) + ' bits of nosync\n')
        row[n] += 1
        col[n] = 1
        nosync = 0
      if framestart == stream_end[n]:
        stream_end[n] = 0
      else: # there is something behind the flush - store it
        print('something behind flush')
        j = 0
        for i in range(framestart+1, stream_end[n]+1):
          j += 1
          stream[n][j] = stream[n][i]
        stream_end[n] = j
      return
    # find nosynch of frames
    # count and skip 1
    for framestart in range(framestart, stream_end[n]+1):
      if stream[n][framestart] == 1:
        nosync += 1
      else:
        break
    else:    # if at the end of the stream store the nosync count to continue next time and return
      if nosync:
        stream[n][1] = -nosync # store number of nosync bytes
        stream_end[n] = 1
      else:
        stream_end[n] = 0
        print('this is unexpected')
        out[n].write('this is unexpected')
      return
    # count 0, skip all but 8
    eight = 0
    for framestart in range(framestart, stream_end[n]+1):
      if stream[n][framestart] == 0:
        eight += 1
      else:
        break
    else:  # at the end of available data
      if eight > 8:
        nosync += eight - 8
        framestart -= 7
      else:
        framestart -= eight - 1
      if nosync:
        stream[n][1] = -nosync # store number of nosyn bytes
        j = 1
      else:
        j = 0
      # leave bytes from posible start of frame
      for i in range(framestart, stream_end[n]+1):
        j += 1
        stream[n][j] = stream[n][i]
      stream_end[n] = j
      return
    if eight < 8:
      nosync += eight
      continue
    else:
      nosync += eight - 8
      framestart -= 8
      eight = 8
      # check for odd bytes within possible frame
      for i in range(framestart+8, stream_end[n]+1):
        if stream[n][i] != 1 and stream[n][i] != 0:
          # found, increase framestart and nosync and break this for
          nosync += i - framestart
          framestart = i
          break
        if i >= framestart + 79:
          break
      if stream[n][i] != 1 and stream[n][i] != 0:
        continue #continue outside loop
      if i < framestart + 79: # not enough bytes to do something usefull, cleanup and return
        if nosync:
          stream[n][1] = -nosync # store number of nosyn bytes
          j = 1
        else:
          j = 0
        # leave bytes from posible start of frame
        for i in range(framestart, stream_end[n]+1):
          j += 1
          stream[n][j] = stream[n][i]
        stream_end[n] = j
        return
      # check for other 9 other frame synch bytes
      for i in range(1, 10):                                          
        if  stream[n][framestart+i*8] != 1: # not a frame
          if i == 1:
            print('unexpected 0 as this byte is already checked to be 1\n')
            out[n].write('unexpected 0 as this byte is already checked to be 1\n')
            framestart += 8
            nosync += 8
          else:
            # right behind last found 1 is first chance for finding 8 x 0
            framestart += (i-1)*8 + 1
            nosync += (i-1)*8 + 1
          break
      else: # end of for, no break so frame found
        if nosync: # end of nosynch found
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(nosync) + ' bits of nosync \n')
          row[n] += 1
          col[n] = 1
          nosync = 0
        for j in range(1, 81):
          frame[j] = stream[n][framestart + j - 1]
        check_v110_frame(frame, n)
        framestart += 80
  # clean out of while loop, all bytes handled, zero bytes left
  stream_end[n] = 0 # just ended at the end of a frame

flush = 2
def rate_adaption(bslen, bytestream, n):
  bits = []
  nbits = 0
  # check validity of bytes
  for i in range(0, bslen):
    if count[n] and bytestream[i] == samebyte[n]:
      count[n] += 1
      continue
    if adaption_rate == 8:
      if (mask_check and bytestream[i] & 0x7F) != 0x7F:
        #flush
        bits.append(flush); nbits += 1
        adapted_find_v110_frames(nbits, bits, n)      
        bits = [] ; nbits = 0
        if count[n]:
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(count[n]) + ' x 0x' + '{:02X}'.format(samebyte[n]) + ' bytes\n');
          row[n] += 1
          col[n] = 1
        samebyte[n] = bytestream[i]
        count[n] = 1
      else:
        bits.append((bytestream[i] >> 7) & 1); nbits += 1
        if count[n]:
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(count[n]) + ' x 0x' + '{:02X}'.format(samebyte[n]) + ' bytes\n');
          row[n] += 1
          col[n] = 1
          count[n] = 0
    elif adaption_rate == 16:
      if (mask_check and bytestream[i] & 0x3F) != 0x3F:
        #flush
        bits.append(flush); nbits += 1
        adapted_find_v110_frames(nbits, bits, n)      
        bits = [] ; nbits = 0
        if count[n]:
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(count[n]) + ' x 0x' + '{:02X}'.format(samebyte[n]) + ' bytes\n');
          row[n] += 1
          col[n] = 1
        samebyte[n] = bytestream[i]
        count[n] = 1         
      else:
        bits.append((bytestream[i] >> 7) & 1); nbits += 1
        bits.append((bytestream[i] >> 6) & 1); nbits += 1
        if count[n]:
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(count[n]) + ' x 0x' + '{:02X}'.format(samebyte[n]) + ' bytes\n');
          row[n] += 1
          col[n] = 1
          count[n] = 0
    elif adaption_rate == 32:
      if (mask_check and bytestream[i] & 0x0F) != 0x0F:
        #flush
        bits.append(flush); nbits += 1
        adapted_find_v110_frames(nbits, bits, n)      
        bits = [] ; nbits = 0
        if count[n]:
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(count[n]) + ' x 0x' + '{:02X}'.format(samebyte[n]) + ' bytes\n');
          row[n] += 1
          col[n] = 1
        samebyte[n] = bytestream[i]
        count[n] = 1         
      else:
        bits.append((bytestream[i] >> 7) & 1); nbits += 1
        bits.append((bytestream[i] >> 6) & 1); nbits += 1
        bits.append((bytestream[i] >> 5) & 1); nbits += 1
        bits.append((bytestream[i] >> 4) & 1); nbits += 1
        if count[n]:
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(count[n]) + ' x 0x' + '{:02X}'.format(samebyte[n]) + ' bytes\n');
          row[n] += 1
          col[n] = 1
          count[n] = 0
    elif adaption_rate == 64:
      bits.append((bytestream[i] >> 7) & 1); nbits += 1
      bits.append((bytestream[i] >> 6) & 1); nbits += 1
      bits.append((bytestream[i] >> 5) & 1); nbits += 1
      bits.append((bytestream[i] >> 4) & 1); nbits += 1
      bits.append((bytestream[i] >> 3) & 1); nbits += 1
      bits.append((bytestream[i] >> 2) & 1); nbits += 1
      bits.append((bytestream[i] >> 1) & 1); nbits += 1
      bits.append(bytestream[i] & 1); nbits += 1
  adapted_find_v110_frames(nbits, bits, n)      

def read_bytes_from_QATS_file(fn):
  global nx
  global nn
  global bytess
  global stream
  global stream_end
  global stream_max_end
  global bitstream
  global bitstream_end
  global bitstream_max_end
  global bytestream
  global bytestream_end
  global bytestream_max_end
  global out
  global row
  global col
  global idrow
  global samebyte
  global count
  global linetime

  f = open(fn + '.slf', 'rb', 99999)
  if read_t_UI32(f) != 4294901808:
    print('magic number not slf version3: 3000 FFFF')
  fileHeaderLen = read_t_UI32(f)
  initial_Timestamp = read_t_UI64(f)
  file_Timestamp = read_t_UI64(f)
  last_Relative_Timestamp = read_t_UI32(f)
  records_in_file = read_t_UI32(f)
  len = read_t_UI16(f)
  for a in range(0,  len):
    read_t_UI8(f)
  extra_length = read_t_UI32(f)
  for a in range(0,  extra_length):
    read_t_UI8(f)
  for a in range(0, records_in_file):
    record_length = read_t_UI16(f)
    source_Id = read_t_UI16(f)
    record_number = read_t_UI32(f)
    relative_timestamp = read_t_UI32(f)
    extra_length = read_t_UI8(f)
    for b in range(0, extra_length):
      read_t_UI8(f)
    n = '{:04X}'.format(source_Id)
    if a == 0:
      bytess = list(range(200))
      for n in out:
        linetime[n] = 0
    if (record_length - 13 - extra_length == 68):
      if not (n in out):
        row[n] = 1
        col[n] = 1
        idrow[n] = 0
        stream[n] = []
        stream[n].append(0)
        stream_end[n] = 0
        stream_max_end[n] = 0
        bitstream[n] = []
        bitstream_end[n] = 0
        bitstream_max_end[n] = 0
        bytestream[n] = []
        bytestream_end[n] = 0
        bytestream_max_end[n] = 0
        out[n] = open(fn + '_' + n + '.txt', 'w')
        samebyte[n] = -1
        count[n] = 0
        if bit_rate%4800 == 0:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t011 E4-E7\tD1-D48                                          \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        elif bit_rate == 2400:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t110 E4-E7\tD1-D24                  \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        elif bit_rate == 1200:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t010 E4-E7\tD1-D12      \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        elif bit_rate == 600:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t100 E4-E7\tD1-D6 \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        if decode_ascii:
          out[n].write('\tASCII')
        if decode_hdlc:
          out[n].write('\tDHCL')
        out[n].write('\n')          
        nx += 1
        nn[nx] = n
        linetime[n] = 0
      for i in range(0, record_length - 13 - extra_length):
        bytess[i] = read_t_UI8(f)
      if (relative_timestamp > linetime[n]): # filter out weird stuff
        if row[n] != idrow[n]:
          out[n].write(str(record_number) + '\t' + str(datetime.datetime.fromtimestamp((file_Timestamp + relative_timestamp)//1000).strftime('%Y-%m-%d %H:%M:%S')) + '.%03d' % ((file_Timestamp + relative_timestamp)%1000) + '\t')
          col[n] = 3
          idrow[n] = row[n]
        rate_adaption(record_length - 13 - extra_length, bytess, n)
        linetime[n] = relative_timestamp
      #else:
        #print(n + '\t' + str(record_number) + '\t' + str(datetime.datetime.fromtimestamp((file_Timestamp + relative_timestamp)//1000).strftime('%Y-%m-%d %H:%M:%S')) + '.%03d' % ((file_Timestamp + relative_timestamp)%1000) + '\t' + str(datetime.datetime.fromtimestamp((file_Timestamp + linetime[n])//1000).strftime('%Y-%m-%d %H:%M:%S')) + '.%03d' % ((file_Timestamp + linetime[n])%1000))			
    else: # still read but do not store
      for a in range(0, record_length - 13 - extra_length):
        read_t_UI8(f)
  f.close()

def read_bytes_from_NetHawk_file(fn):
  global nx
  global nx
  global nn
  global bytess
  global stream
  global stream_end
  global stream_max_end
  global bitstream
  global bitstream_end
  global bitstream_max_end
  global bytestream
  global bytestream_end
  global bytestream_max_end
  global out
  global row
  global col
  global idrow
  global samebyte
  global count
  f = open(fn + '.txt', 'r')
  while True:
    for line in f:
      if 'PCM ' in line:
        break
    else: # end of file
      f.close()
      return
    d1, n, d2, d3, ctype, logid, date, time = line.split()
    if ctype == "Type:tfomessage":
      if (line[0] == 'P'):
        n += "_R1"
      else:
        n += "_R2"
      if not (n in out):
        row[n] = 1
        col[n] = 1
        idrow[n] = 0
        stream[n] = []
        stream[n].append(0)
        stream_end[n] = 0
        stream_max_end[n] = 0
        bitstream[n] = []
        bitstream_end[n] = 0
        bitstream_max_end[n] = 0
        bytestream[n] = []
        bytestream_end[n] = 0
        bytestream_max_end[n] = 0
        out[n] = open(fn + '_' + n + '.txt', 'w')
        samebyte[n] = -1
        count[n] = 0
        if bit_rate%4800 == 0:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t011 E4-E7\tD1-D48                                          \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        elif bit_rate == 2400:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t110 E4-E7\tD1-D24                  \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        elif bit_rate == 1200:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t010 E4-E7\tD1-D12      \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        elif bit_rate == 600:
          out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t100 E4-E7\tD1-D6 \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t')
        if decode_ascii:
          out[n].write('\tASCII')
        if decode_hdlc:
          out[n].write('\tDHCL')
        out[n].write('\n')          
        nx += 1
        nn[nx] = n
      if row[n] != idrow[n]:
        out[n].write(str(logid) + '\t' + date[5:] + " " + time + '\t')
        col[n] = 3
        idrow[n] = row[n]
      speech = f.readline()
      if 'SPEECH' in speech:
        x = int(speech.split('(')[1].split()[0])
        if x % 32:
          lines = x//32+1
        else:
          lines = x//32
        bytess = []
        for i in range(0, lines):
          line = f.readline()
          bytess.extend(map(int, bytes.fromhex(line.strip().replace(' ',''))))
        rate_adaption(len(bytess), bytess, n)


        
""" 
stream = {}
stream_end = {}
stream_max_end = {}
bitstream = {}
bitstream_end = {}
bitstream_max_end = {}
bytestream = {}
bytestream_end = {}
bytestream_max_end = {}
out = {}
row = {}
col = {}
idrow = {}
samebyte = {}
count = {}
bit_rate = 9600#bit/s
adaption_rate = 16#kbit/s
read_bytes_from_NetHawk_file('9600')
"""
if len(argv) < 2 or (len(argv) == 2 and argv[1][0] == "-"):
  exit('Usage: %s [-<baudrate>] <filename> [<filename> [<filename> [..]]]' % argv[0])
  print ('No baudrate specified. Using default 4800.')
bit_rate = 4800#baud
adaption_rate = 8#kb/s
stitch = False
decode_ascii = False
decode_hdlc = False
mask_check = True

for arg in range(1, len(argv)):
  if not stitch: 
    # globals
    nx = 0
    nn = {}
    bytess = []
    linetime = {}
    stream = {}
    stream_end = {}
    stream_max_end = {}
    bitstream = {}
    bitstream_end = {}
    bitstream_max_end = {}
    bytestream = {}
    bytestream_end = {}
    bytestream_max_end = {}
    out = {}
    row = {}
    col = {}
    idrow = {}
    samebyte = {}
    count = {}
  if  argv[arg][0] == "-":
    if argv[arg] == '-stitch': # stitch the input files togeter and handle as 1 continues input 
      stitch = True
    elif argv[arg] == '-ascii': # decode payload as ASCII
      decode_ascii = True
    elif argv[arg] == '-hdlc': # decode payload as HDLC
      decode_hdlc = True
    elif argv[arg] == '-nomaskcheck': # ignore masked out bits in rate adaptation
      mask_check = False      
    elif -int(argv[arg]) in [50, 75, 110, 150, 200, 300, 600]:
      bit_rate = 600
      adaption_rate = 8
      #print("Using 600 V.110 bit rate within 8k to 64k bit rate adaption.")
    elif -int(argv[arg]) == 1200:
      bit_rate = 1200
      adaption_rate = 8
      #print("Using 1200 V.110 bit rate within 8k to 64k bit rate adaption.")
    elif -int(argv[arg]) == 2400:
      bit_rate = 2400
      adaption_rate = 8
      #print("Using 3400 V.110 bit rate within 8k to 64k bit rate adaption.")
    elif -int(argv[arg]) in [3600, 4800]:
      bit_rate = 4800
      adaption_rate = 8
      #print("Using n x 4800 V.110 bit rate within 8k to 64k bit rate adaption.")
    elif -int(argv[arg]) in [7200, 9600]:
      bit_rate = 9600
      adaption_rate = 16
      #print("Using n x 4800 V.110 bit rate within 16k to 64k bit rate adaption.")
    elif -int(argv[arg]) in [12000, 14400, 19200]:
      bit_rate = 19200
      adaption_rate = 32
      #print("Using n x 4800 V.110 bit rate within 32k to 64k bit rate adaption.")
    elif -int(argv[arg]) in [24000, 28800, 38400]:
      bit_rate = 38400
      adaption_rate = 64
      #print("Using n x 4800 V.110 bit rate within 64k to 64k bit rate adaption.")
    else:
      print ('%s specifies a unknown option. Using default 4800.' % argv[arg])
  else:         
    filename, fileext = splitext(argv[arg])
    if fileext == '.slf':
      print('decoding %s as QATS protocol analyzer file' % argv[arg])
      read_bytes_from_QATS_file(filename)
    elif fileext == '.txt':
      print('decoding %s as NetHawk protocol analyzer text output' % argv[arg]) 
      read_bytes_from_NetHawk_file(filename)
    else:
      print ('%s has a unknown file extention. Skipping.' % argv[arg])
  if not stitch: 
    # flush
    print('end\n');
    bytess = list(range(80))
    for i in range(0, 80):
      bytess[i] = 0
    for i in range(1, nx+1):
      rate_adaption(80, bytess, nn[i])
      out[nn[i]].close()
if stitch: 
  # flush
  print('end\n');
  bytess = list(range(80))
  for i in range(0, 80):
    bytess[i] = 0
  for i in range(1, nx+1):
    rate_adaption(80, bytess, nn[i])
    out[nn[i]].close()

  
