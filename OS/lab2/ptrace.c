#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/ptrace.h> 
#include <sys/user.h> 
#include "syscallnames.h"

const char* scnum_to_name(int num) {
	return syscallnames[num];
}

int main(int argc, char* argv[]) {
	int status; 
 	struct user_regs_struct regs;
 	int counter = 0, is_entry = 1;
 	pid_t pid = fork();
	if (pid < 0) {
     		perror("Couldn't fork");
     		exit(1);
	}
   	else if (pid == 0) { /* in the child process */
     		ptrace(PTRACE_TRACEME, 0, NULL, NULL);
     		execvp(argv[1], argv+1);
   	}
	else { /* in the parent process */
     		wait(&status);
     		while(status == 1407) {
       			ptrace(PTRACE_GETREGS, pid, NULL, &regs);
       			if (is_entry) {
         			printf("System Call %s (%lld) called with %lld, %lld, %lld\n", 
					scnum_to_name(regs.orig_rax), regs.orig_rax, regs.rdi, regs.rsi, regs.rdx);
         			is_entry = 0;
         			counter++;
       			}
       			else
         			is_entry = 1; 
     			ptrace(PTRACE_SYSCALL, pid, NULL, NULL); 
     			wait(&status); 
     		}
   	}
   	printf("Total Number of System Calls = %d\n", counter);
   	return 0; 
}
