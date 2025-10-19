#ifndef ROVACA_HC_OBJECT_POOL_H_
#define ROVACA_HC_OBJECT_POOL_H_
#include <list>
#include <mutex>

/*!
 * @brief
 * @tparam Type 必须是指针
 */
template <typename Type>
class ObjectPool
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::mutex objs_lock_;
    std::list<Type> objs_;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ObjectPool()
        : objs_lock_()
        , objs_()
    {}

    ~ObjectPool()
    {
        for (Type t : objs_) {
            delete t;
        }
    }

    /*!
     * @brief 有元素返回链表头, 无元素返回nullptr, 由调用者自己new, 调用者自己new的也要返回到此pool, pool 最大容量 = 消费者数量
     * @return
     */
    Type pop()
    {
        Type t = nullptr;
        {
            std::lock_guard<std::mutex> lock(objs_lock_);
            if (!objs_.empty()) {
                t = objs_.front();
                objs_.pop_front();
            }
        }
        return t;
    }

    void push(Type t)
    {
        std::lock_guard<std::mutex> lock(objs_lock_);
        objs_.push_back(t);
    }
    std::list<Type>& obj_list() { return objs_; }
};

#endif  // ROVACA_HC_OBJECT_POOL_H_
