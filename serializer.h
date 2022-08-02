#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <filesystem>

struct ISerializer
{
    virtual ~ISerializer()
    {

    }

    template <class T>
    void serializeArraySize(std::vector<T>& p)
    {
        size_t size = p.size();
        serializePreallocatedMem(&size, sizeof(size));
        p.resize(size);
    }

    template <class T>
    void serializeStdArray(std::vector<T>& p)
    {
        serializeArraySize(p);
        if (p.size() > 0)
        {
            serializePreallocatedMem(&p[0], sizeof(p[0]) * (NvU32)p.size());
        }
    }

    template <class T>
    void serializeArrayOfPointers(std::vector<T*>& p)
    {
        // serialize array size
        serializeArraySize(p);

        // serialize each pointer in the array
        for (NvU32 u = 0; u < p.size(); ++u)
        {
            NvU32 exists = (p[u] != nullptr);
            serializePreallocatedMem(&exists, sizeof(exists));
            if (!exists)
            {
                continue;
            }
            if (m_isReading)
            {
                p[u] = new T();
            }
            p[u]->serialize(*this);
        }
    }

    template <class T>
    void serializeSimpleType(T& value)
    {
        serializePreallocatedMem(&value, sizeof(value));
    }
    virtual void serializePreallocatedMem(void* pMem, NvU32 memSizeInBytes) = 0;

protected:
    bool m_isReading = false;
    std::wstring makeFullPath(const std::filesystem::path& path)
    {
        if (path.is_absolute())
            return path;
        return std::wstring(L"c:\\atomNets\\") + path.c_str();
    }
};

struct MyWriter : public ISerializer
{
    MyWriter(const std::filesystem::path& path)
    {
        std::wstring sFullPath = makeFullPath(path);
        _wfopen_s(&m_fp, sFullPath.c_str(), L"wb");
        nvAssert(m_fp != nullptr);
    }
    ~MyWriter()
    {
        fclose(m_fp);
    }
    virtual void serializePreallocatedMem(void* pMem, NvU32 memSizeInBytes) override
    {
        nvAssert(m_fp != nullptr && memSizeInBytes > 0);
        fwrite(pMem, 1, memSizeInBytes, m_fp);
    }
private:
    FILE* m_fp = nullptr;
};

struct MyReader : public ISerializer
{
    MyReader(const std::filesystem::path& path)
    {
        std::wstring sFullPath = makeFullPath(path);
        _wfopen_s(&m_fp, sFullPath.c_str(), L"rb");
        nvAssert(m_fp != nullptr);
        m_isReading = true;
    }
    ~MyReader()
    {
        fclose(m_fp);
    }
    virtual void serializePreallocatedMem(void* pMem, NvU32 memSizeInBytes) override
    {
        nvAssert(m_fp != nullptr);
        fread(pMem, 1, memSizeInBytes, m_fp);
    }
private:
    FILE* m_fp = nullptr;
};