.intel_syntax noprefix

.text
.global bn_muls_low

bn_muls_low:
    push   r12
    mov    rax, rcx
    shr    rax, 63
    xor    rax, rdx     # sign = sa ^ sd;
    mov    rdx, rcx
    neg    rcx
    cmovns rdx, rcx     # rdx = (digit < 0 ? -digit : digit);
    xor    rcx, rcx     # clear flags, create zero

    mulx   r9, r8, [rsi]

    mulx   r10, r11, [rsi + 8]
    adcx   r9, r11

    mulx   r11, r12, [rsi + 16]
    adcx   r10, r12

    mulx   r12, rsi, [rsi + 24]
    adcx   r11, rsi
    adcx   r12, rcx

    neg    rax
    xor    r8, rax
    xor    r9, rax
    xor    r10, rax
    xor    r11, rax
    xor    r12, rax

    neg    rax
    add    r8, rax
    adcx   r9, rcx
    adcx   r10, rcx
    adcx   r11, rcx
    adcx   r12, rcx
    mov    [rdi], r8
    mov    [rdi + 8], r9
    mov    [rdi + 16],r10
    mov    [rdi + 24],r11
    mov    rax, r12

    pop    r12
	ret

.global bn_addn_low
bn_addn_low:
	xor    rax, rax
	mov    r8, [rdx]
	add    r8, [rsi]
	mov    [rdi], r8
	mov    r9, [rdx + 8]
	adc    r9, [rsi + 8]
	mov    [rdi+8], r9
	mov    r10, [rdx + 16]
	adc    r10, [rsi + 16]
	mov    [rdi + 16], r10
	mov    r11, [rdx + 24]
	adc    r11, [rsi + 24]
	mov    [rdi + 24], r11
	mov    rax, [rdx + 32]
	adc    rax, [rsi + 32]
	mov    [rdi + 32], rax
	ret