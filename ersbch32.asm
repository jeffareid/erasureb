	title	ersbch32
;	page	80,132
;-----------------------------------------------------------------------;
;	ersbch32.asm	fast GF(2^8) stuff				;
;		    Copyright(c) 2020, Jeff Reid			;
;		    2020JUN27 21:15					;
;-----------------------------------------------------------------------;
nrow	equ	20			;max # rows
mx	struct				;matrix structure
p	dq	?			;ptr to bfr
m	dq	nrow dup (?)		;ptr to matrix
r	dq	?			;# rows
c	dq	?			;# columns
mx	ends
;-----------------------------------------------------------------------;
;	DATA								;
;-----------------------------------------------------------------------;
	.data
	align 16
mask1	dq	00f0f0f0f0f0f0f0fh
	dq	00f0f0f0f0f0f0f0fh
mask2	dq	0f0f0f0f0f0f0f0f0h
	dq	0f0f0f0f0f0f0f0f0h
	.data?
	align 16
mpytbl	db	(256*2*16) dup (?)	;mpy tables

	.code
;	c++ name mangled
	public	?InitGFa@@YAXPEAE@Z
	align 16
;	InitGFa		rcx = BYTE *abMpy
?InitGFa@@YAXPEAE@Z proc
	lea	rdx,mpytbl
;	copy	mpy by 0x01 
	lea	r9,[rcx+65536]
init0:	lea	r8,[rcx+16]
init1:	mov	al,[rcx]
	mov	[rdx],al
	inc	rcx
	inc	rdx
	cmp	rcx,r8
	jne	short init1
	add	rcx,256-16
	add	rdx,16
	cmp	rcx,r9
	jb	short init0
	sub	rcx,65536
;	copy	mpy by 0x10
	lea	rdx,mpytbl+16
init2:	lea	r8,[rcx+256]
init3:	mov	al,[rcx]
	mov	[rdx],al
	add	rcx,16
	inc	rdx
	cmp	rcx,r8
	jne	short init3
	add	rdx,16
	cmp	rcx,r9
	jb	short init2
	ret
?InitGFa@@YAXPEAE@Z endp

;	c++ name mangled
	public	?MatrixMpya@@YAXAEAVMATRIX@@00@Z
	align 16
;	MatrixMpya	rcx = MATRIX &mDst
;			rdx = MATRIX &mSrc0
;			r8  = MATRIX &mSrc1
?MatrixMpya@@YAXAEAVMATRIX@@00@Z proc
	push	rbp
	mov	rbp,rsp
	and	rsp,0fffffffffffffff0h
	push	r12
	push	r13
	push	r14
	push	r15
	sub	rsp,96
	movdqa	[rsp+16*0],xmm6
	movdqa	[rsp+16*1],xmm7
	movdqa	[rsp+16*2],xmm8
	movdqa	[rsp+16*3],xmm9
	movdqa	[rsp+16*4],xmm10
	movdqa	[rsp+16*5],xmm11
	movdqa	xmm10,xmmword ptr mask1
	movdqa	xmm11,xmmword ptr mask2
	lea	r9,mpytbl
	xor	r15,r15			;r15 = r = dst row index
mpy0:	mov	r12,[rcx+r15*8].mx.m	;r12 =	dst row ptr
	mov	r11,[rdx+r15*8].mx.m	;r11 = src0 row ptr
	xor	r14,r14			;r14 = m = inner dst col = src0 row index
mpy1:	mov	r10,[r8 +r14*8].mx.m	;r10 = src1 row ptr
	movzx	rax,byte ptr[r11+r14]	;rax = b = multiplier
	shl	rax,5			;rax = index into mpytbl
	movdqa	xmm8,xmmword ptr [rax+r9+ 0] ;xmm8 = mpy lo nibble by b tbl
	movdqa	xmm9,xmmword ptr [rax+r9+16] ;xmm9 = mpy hi nibble by b tbl
	xor	r13,r13			;r13 = c = outer dst col = src1 col index
mpy2:	movdqa	xmm4,xmmword ptr [   r10+r13] ;get 32 bytes
	movdqa	xmm5,xmmword ptr [16+r10+r13]
	movdqa	xmm6,xmm4		;copy data to split into hi + lo nibbles
	movdqa	xmm7,xmm5
	pand	xmm4,xmm10		;mask data into nibbles
	pand	xmm5,xmm10
	pand	xmm6,xmm11
	pand	xmm7,xmm11
	psrlq	xmm6,4			;shift hi nibbles to lo nibbles
	psrlq	xmm7,4
	movdqa	xmm0,xmm8		;copy hi+lo table data
	movdqa	xmm1,xmm8
	movdqa	xmm2,xmm9
	movdqa	xmm3,xmm9
	pshufb	xmm0,xmm4		;mpy nibbles by constant via pshufb
	pshufb	xmm1,xmm5
	pshufb	xmm2,xmm6
	pshufb	xmm3,xmm7
	pxor	xmm0,xmm2		;xor the two sets of hi+lo products
	pxor	xmm1,xmm3
	or	r14,r14			;br if first time (just store)
	jz	short mpy3		; else xor
	pxor	xmm0,xmmword ptr [   r12+r13]
	pxor	xmm1,xmmword ptr [16+r12+r13]
mpy3:	movdqa	xmmword ptr [	r12+r13],xmm0
	movdqa	xmmword ptr [16+r12+r13],xmm1
	add	r13,32
	cmp	r13,[ r8].mx.c
	jb	mpy2
	inc	r14
	cmp	r14,[rdx].mx.c
	jb	mpy1
	inc	r15
	cmp	r15,[rdx].mx.r
	jb	mpy0
	movdqa	xmm6, [rsp+16*0]
	movdqa	xmm7, [rsp+16*1]
	movdqa	xmm8, [rsp+16*2]
	movdqa	xmm9, [rsp+16*3]
	movdqa	xmm10,[rsp+16*4]
	movdqa	xmm11,[rsp+16*5]
	add	rsp,96
	pop	r15
	pop	r14
	pop	r13
	pop	r12
	mov	rsp,rbp
	pop	rbp
	ret
?MatrixMpya@@YAXAEAVMATRIX@@00@Z endp

;	c++ name mangled
	public	?MatrixXora@@YAXAEAVMATRIX@@_K@Z
	align 16
;	MatrixXora	rcx = MATRIX &mDst
;			rdx = row, != 0
?MatrixXora@@YAXAEAVMATRIX@@_K@Z proc
	push	rdi
	push	rsi
	mov	r11,[rcx].mx.r		;r11  = # rows
	mov	r10,[rcx].mx.c		;r10  = # cols
	mov	r9,32			;r9   = # cols per read/write
	mov	rdi,[rcx+rdx*8].mx.m	;rdi = ptr to dst row
	mov	rsi,[rcx].mx.m		;rsi = ptr to row 0
	mov	r8,1			;r8 = row index
	cmp	rsi,rdi			;br if dst row != row 0
	jne	short xor0
	mov	rsi,[rcx+r8*8].mx.m	;rsi = ptr to row 1
	mov	r8,2			;r8 = row index 2
xor0:	lea	rax,[rsi+r10]		;rax = end 1st row
xor1:	movdqa	xmm0,xmmword ptr [rsi]	;copy src row to dst row
	movdqa	xmm1,xmmword ptr [16+rsi]
	add	rsi,r9
	movdqa	xmmword ptr [rdi],xmm0
	movdqa	xmmword ptr [rdi+16],xmm1
	add	rdi,r9
	cmp	rsi,rax
	jb	xor1
xor2:	cmp	r8,rdx			;br if dst row
	je	short xor4
	mov	rsi,[rcx+r8*8].mx.m	;rsi = ptr to src row
	mov	rdi,[rcx+rdx*8].mx.m	;rdi = ptr to dst row
	lea	rax,[rsi+r10]		;rax = end src row
xor3:	movdqa	xmm0,xmmword ptr [rsi]	;xor src to dst
	movdqa	xmm1,xmmword ptr [rsi+16]
	add	rsi,r9
	pxor	xmm0,xmmword ptr [rdi]
	pxor	xmm1,xmmword ptr [rdi+16]
	movdqa	xmmword ptr [rdi],xmm0
	movdqa	xmmword ptr [rdi+16],xmm1
	add	rdi,r9
	cmp	rsi,rax
	jb	xor3
xor4:	inc	r8
	cmp	r8,r11
	jb	xor2
	pop	rsi
	pop	rdi
	ret
?MatrixXora@@YAXAEAVMATRIX@@_K@Z endp

	end

