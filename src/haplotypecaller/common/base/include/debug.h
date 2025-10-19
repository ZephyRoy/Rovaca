#ifndef __DEBUG_H__
#define __DEBUG_H__

typedef void (*p_rovaca_log_func)(int, const char*, const char*, const char*);
extern p_rovaca_log_func rovaca_assemble_log_trace;
extern p_rovaca_log_func rovaca_assemble_log_debug;
extern p_rovaca_log_func rovaca_assemble_log_info;
extern p_rovaca_log_func rovaca_assemble_log_warn;
extern p_rovaca_log_func rovaca_assemble_log_error;
extern p_rovaca_log_func rovaca_assemble_log_critical;

#define statics_print_trace(msg) rovaca_assemble_log_trace(__LINE__, __FILE__, __func__, msg);
#define statics_print_debug(msg) rovaca_assemble_log_debug(__LINE__, __FILE__, __func__, msg);
#define statics_print_info(msg)  rovaca_assemble_log_info(__LINE__, __FILE__, __func__, msg);
#define statics_print_warn(msg)  rovaca_assemble_log_warn(__LINE__, __FILE__, __func__, msg);
#define statics_print_log(msg)   rovaca_assemble_log_error(__LINE__, __FILE__, __func__, msg);
#define statics_print_error(msg) rovaca_assemble_log_critical(__LINE__, __FILE__, __func__, msg);

#endif  // !__DEBUG_H__