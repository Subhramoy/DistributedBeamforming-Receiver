import socket
import struct
import sys
import threading
import logging
import time



## This class should not be in use,
#  customized version of this is defined in CSI_fb_adapter module
class UDPServer(threading.Thread):
    def __init__(self, udp_packet_handler):
        print("UDP server generated.")
        self.udp_packet_handler = udp_packet_handler
        threading.Thread.__init__(self)

    def run(self):
        multicast_group = '224.3.29.71'
        server_address = ('', 10000)

        # Create the socket
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

        # Bind to the server address
        sock.bind(server_address)

        # Tell the operating system to add the socket to the multicast group
        # on all interfaces.
        group = socket.inet_aton(multicast_group)
        mreq = struct.pack('4sL', group, socket.INADDR_ANY)
        sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

        # Receive/respond loop
        while True:
            print >> sys.stderr, '\nwaiting to receive message'
            data, address = sock.recvfrom(1024)
            print(data)
            self.csi.set_CSI(data)

            print >> sys.stderr, 'received %s bytes from %s' % (len(data), address)
            print >> sys.stderr, data

            print >> sys.stderr, 'sending acknowledgement to', address
            sock.sendto('ack', address)

