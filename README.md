# V.110-decoder
Python3 decoder script for raw 4800 baud transparant data PCM data to V.110 frames and content

v.110.py reads files and interprets the data in each file as V.110 data
By file extenttion will tell the type of file. At the moment .tt is interpreted as NetHawk text output and .slf as QATS binairy output

The baudrate can be set Usage: python v.110.py [-\<baudrate\>] \<filename\> [\<filename\>][..]

The next functions are steps refine the data and then call the function below

read_bytes_from_QATS_file(filename) or read_bytes_from_NetHawk_file(filename)
  Reads parts from a file and creates the bytestream to be analysed
  If the cursor is at the beginning of the line a timestamp is written to output

rate_adaption()
  Perfroms rate adpation by filtering the relevant bits from the bytestream
  Byte values which are not used for the specified rates are counted and written to output

adapted_find_v110_frames()
  Finds V.110 frames in the bit stream
  Writes the number of bytes that don't belong to any frame to the output

check_v110_frame()
  Gets the part of the embedded bitstream from the V.110 frames 
  Writes all the V.110 frame bits to the output
  
check_bitstream()
  Turns the embedded data bitstream into the bytes that are send via the V.110 protocol (is not passed on)
  Writes bytes and breaks signals
  Writes bytes converted to ASCII (It would be cleaner if this part is moved to its own function)
  
To add a new file format a new function read_bytes_from_XXX_file has to e written and called from the main function at the bottom.
To add a new decoding of the send bytes write a new function to be called by check_bitstream(),  remove the ASCII part at the bottom of check_bitstream()


