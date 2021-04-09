import os

def write_input(npath, replicas):
	with open("rs.run", "w") as rsfile, open("remd.groupfile","w") as groupfile:
		for i in range(1,replicas+1):
			groupfile.write("-O -rem 3 -l "+npath+"/logfile."+str(i)+" -remlog "+npath+"/rem.log -i remd.in -o "+npath+"/remd.out."+str(i)+" -c "+npath+"/heated.rst -r "+npath+"/remd.rst."+str(i)+" -x "+npath+"/remd."+str(i)+".nc -inf "+npath+"/remd.mdinfo."+str(i)+" -p "+npath+"/system_"+str(i)+".top\n")
		rsfile.write("mpirun --oversubscribe -np "+str(replicas)+" pmemd.cuda.MPI -ng "+str(replicas)+" -groupfile "+npath+"/remd.groupfile\n")


###################
# MAIN
###################

replicas = 4
npath = '.'
write_input(npath, replicas)




































