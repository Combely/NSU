#include <stdio.h>
#include <sys/ptrace.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/reg.h>  

void do_child() {
        ptrace(PTRACE_TRACEME, 0, 0, 0); 
        execl("/bin/echo", "echo", "Hello world!",  NULL);
	perror("execl");
}


void do_parent(pid_t child_pid) {
        int wait_val, counter = 0;
	wait(&wait_val); 
	while (wait_val == 1407 ) {
                counter++;
                if (ptrace(PTRACE_SINGLESTEP, child_pid, 0, 0) != 0)
                	perror("ptrace");
        	wait(&wait_val);
	}
	printf("Number of machine instructions : %d\n", counter);
}


int main() {   
	pid_t pid = fork();
	if (pid < 0) {
		perror("fork");
	}	
	else if (pid == 0) {
    		do_child();
	}
    	else {
		do_parent(pid);
    	}
    	return 0;
}
