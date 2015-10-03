#!/usr/bin/env python
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


def check_bitstream(b, n):
  if bitstream_end[n] + 48 > bitstream_max_end[n]:
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
    

def find_v110_frames(bslen, bytestream, n):
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
    count = 0
    framestart = 2
  else:
    nosync = 0
    count = stream[n][1] >> 8
    framestart = 1
    stream[n][1] = stream[n][1] & 255
    
  while framestart <= stream_end[n]:
    # check stream for x times of same bytes not FF or 7F
    while (stream[n][framestart] != 255 and stream[n][framestart] != 127):
      if nosync:
        for col[n] in range(col[n], 3):
          out[n].write('\t')
        out[n].write(str(nosync) + ' bytes of nosync\n')
        row[n] += 1
        col[n] = 1
        nosync = 0
      samebyte = stream[n][framestart]
      for i in range(framestart+1, stream_end[n]+1):
        if stream[n][i] == samebyte:
          count += 1
        else:
          break
      else:  # run out of bytes, return and continue counting next run
        stream[n][1] = (count << 8) + samebyte
        stream_end[n] = 1
        return
      framestart = i
      # output so far
      for col[n] in range(col[n], 3):
        out[n].write('\t')
      out[n].write(str(count+1) + ' x 0x' + '{:02X}'.format(samebyte) + ' bytes\n');
      row[n] += 1
      col[n] = 1
      count = 0

    # find nosynch of frames

    # count and skip 0xFF
    for framestart in range(framestart, stream_end[n]+1):
      if stream[n][framestart] == 255:
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
      return

    # count 0x7F, skip all but 8
    eight = 0
    for framestart in range(framestart, stream_end[n]+1):
      if stream[n][framestart] == 127:
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
        if stream[n][i] != 255 and stream[n][i] != 127:
          # found, increase framestart and nosync and continue
          nosync += i - framestart
          framestart = i
          continue
        if i >= framestart + 79:
          break
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
        if  stream[n][framestart+i*8] != 255: # not a frame
          if i == 1:
            out[n].write('unexpected 0 as this byte is already checked to be 1\n')
            framestart += (i-1)*8 + 1 # right behind last found FF is first chance for finding 8 x 7F
            nosync += (i-1)*8 + 1
          else:
            framestart += (i-1)*8 + 1 # right behind last found FF is first chance for finding 8 x 7F
            nosync += (i-1)*8 + 1
          break
      else: # end of for, no break so frame found
        if nosync: # end of nosynch found
          for col[n] in range(col[n], 3):
            out[n].write('\t')
          out[n].write(str(nosync) + ' bytes of nosync \n')
          row[n] += 1
          col[n] = 1
          nosync = 0
        for j in range(1, 81):
          frame[j] = stream[n][framestart + j - 1] >> 7
        check_v110_frame(frame, n)
        framestart += 80
  stream_end[n] = 0 # just ended at the end of a frame



def read_bytes_from_QATS_file(fn):
  nx = 0
  nn = {}
  bytess = []
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
    record_lengt = read_t_UI16(f)
    source_Id = read_t_UI16(f)
    record_number = read_t_UI32(f)
    relative_timestamp = read_t_UI32(f)
    extra_length = read_t_UI8(f)
    for b in range(0, extra_length):
      read_t_UI8(f)
    if record_lengt - 13 - extra_length == 68:
      n = '{:04X}'.format(source_Id)
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
        out[n] = open(fn + '_' + n + '.txt', 'w')
        out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t011 E4-E7\tD1-D48                                          \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t\tascii\n')
        nx += 1
        nn[nx] = n
      if row[n] != idrow[n]:
        out[n].write(str(record_number) + '\t' + str(datetime.datetime.fromtimestamp((file_Timestamp + relative_timestamp)//1000).strftime('%Y-%m-%d %H:%M:%S')) + '.%03d' % ((file_Timestamp + relative_timestamp)%1000) + '\t')
        col[n] = 3
        idrow[n] = row[n]
      for i in range(0, 68):
        if a == 0:
          bytess.append(read_t_UI8(f))
        else:  
          bytess[i] = read_t_UI8(f)
      find_v110_frames(68, bytess, n)
    else: # still read but not store
      for a in range(0, record_lengt - 13 - extra_length):
        read_t_UI8(f)
  print('\nend\n');
  for i in range(0, 68):
    bytess[i] = 0
  for i in range(1, nx+1):
    find_v110_frames(68, bytess, nn[i])
    out[nn[i]].close()
  f.close()

def read_bytes_from_NetHawk_file(fn):
  nx = 0
  nn = {}
  bytess = []
  
  f = open(fn + '.txt', 'r')
  while True:
    for line in f:
      if 'PCM ' in line:
        break
    else: # end of file
      print("\nend\n")
      # flush
      for i in range(0, 80):
        bytess[i] = 0
      for i in range(1, nx+1):
        find_v110_frames(68, bytess, nn[i])
        out[nn[i]].close()
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
        out[n] = open(fn + '_' + n + '.txt', 'w')
        out[n].write('log id\ttime\tdata kind\tSA 1.3.6.8\tSB 4.9\tX 2.7\t011 E4-E7\tD1-D48                                          \tbyte1\tbyte2\tbyte3\tbyte4\tbyte5\t\t\tascii\n')
        nx += 1
        nn[nx] = n
      if row[n] != idrow[n]:
        out[n].write(str(logid) + '\t' + date[5:] + " " + time + '\t')
        if logid == 'Id:18974':
          print(logid)
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
        find_v110_frames(len(bytess), bytess, n)
"""
stream = {}
stream_end = {}
stream_max_end = {}
bitstream = {}
bitstream_end = {}
bitstream_max_end = {}
out = {}
row = {}
col = {}
idrow = {}
read_bytes_from_NetHawk_file('NetHawk_trace')
"""
if len(argv) < 2:
    exit('Usage: %s <filename> [<filename> [<filename> [..]]]' % argv[0])
for arg in range(1, len(argv)):
  # globals
  stream = {}
  stream_end = {}
  stream_max_end = {}
  bitstream = {}
  bitstream_end = {}
  bitstream_max_end = {}
  out = {}
  row = {}
  col = {}
  idrow = {}
  filename, fileext = splitext(argv[arg])
  if fileext == '.slf':
    print('decoding %s as QATS protocol analyzer file' % argv[arg])
    read_bytes_from_QATS_file(filename)
  elif fileext == '.txt':
    print('decoding %s as NetHawk protocol analyzer text output' % argv[arg]) 
    read_bytes_from_NetHawk_file(filename)
  else:
    print ('%s has a unknown file extention. Skipping.' % argv[arg])
    
